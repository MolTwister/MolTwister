//
// Copyright (C) 2021 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include "SimulationBox.h"
#include "../Integrators/Constants.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

CSimulationBox::CSimulationBox(CMolTwisterState* state, FILE* stdOut, SMolDynConfigStruct config)
    : velVerlet_(state, stdOut, config.cutoffRadius_, config.neighListShell_, config.cutoffForce_)
{
    state_ = state;
    molDynConfig_ = config;

    velVerlet_.setNonBondedScaleFactors((float)molDynConfig_.scale12Interactions_,
                                       (float)molDynConfig_.scale13Interactions_,
                                       (float)molDynConfig_.scale14Interactions_,
                                       (float)molDynConfig_.scaleAbove14BondedInteractions_);
    velVerlet_.setVerboseOutput(molDynConfig_.verboseOutput_);
    NH_T_.setVerboseOutput(molDynConfig_.verboseOutput_);
    NH_P_.setVerboseOutput(molDynConfig_.verboseOutput_);

    N_ = 100;
    C3DRect pbc = state->view3D_->calcPBC();
    LmaxX_ = pbc.getWidthX(); // [AA]
    LmaxY_ = pbc.getWidthY(); // [AA]
    LmaxZ_ = pbc.getWidthZ(); // [AA]
    dt_ = 1.0 / Conv_t;   // [fs]
    dim_ = dim3D;
    resizeArrays();

    initSystem();
}

void CSimulationBox::initSystem()
{
    NH_T_.M_ = molDynConfig_.temperatureNHChainLength_;
    NH_P_.M_ = molDynConfig_.pressureNHChainLength_;

    NH_T_.n_ = molDynConfig_.temperatureRESPASteps_;
    NH_P_.n_ = molDynConfig_.pressureRESPASteps_;
    N_ = (int)state_->atoms_.size();
    NH_T_.T_ = NH_P_.T_ = molDynConfig_.temperature_ / Conv_T;
    dt_ = molDynConfig_.timeStep_ / Conv_t;
    NH_T_.tau_ = molDynConfig_.temperatureRelaxationTime_ / Conv_t;
    NH_P_.tau_ = molDynConfig_.pressureRelaxationTime_ / Conv_t;
    dim_ = dim3D;
    resizeArrays();
    NH_T_.setRandNHPos();
    NH_P_.setRandNHPos();
    setRandMom();
    NH_T_.setRandNHMom();
    NH_P_.setRandNHMom();

    velVerlet_.W_ = NH_P_.T_*NH_P_.tau_*NH_P_.tau_;
    velVerlet_.setRandMom(NH_P_.tau_);
    velVerlet_.P_ = molDynConfig_.pressure_ / Conv_press;
    velVerlet_.V0_ = LmaxX_*LmaxY_*LmaxZ_;

    copyPosFromState();
}

void CSimulationBox::resizeArrays()
{
    NH_T_.prepareArrays(N_, dim_);
    NH_P_.prepareArrays(1, 1);
    particles_.resize(N_);
}

void CSimulationBox::copyPosFromState()
{
    for(int i=0; i<N_; i++)
    {
        if(state_->atoms_[i]->r_.size() == 0) continue;
        particles_[i].m_ = state_->atoms_[i]->m_;
        particles_[i].r_ = state_->atoms_[i]->r_[0];
    }
}

void CSimulationBox::setRandMom()
{
    double  a = 2.0 / double(RAND_MAX);
    
    // Randomize p between 0 and 1 and prevent p=0
    for(int k=0; k<N_; k++)
    {
        particles_[k].p_.x_ = (a*double(rand()) - 1.0);
        particles_[k].p_.y_ = (a*double(rand()) - 1.0);
        particles_[k].p_.z_ = (a*double(rand()) - 1.0);
        
        if(particles_[k].p_.x_ == 0.0) particles_[k].p_.x_ = 1.0;
        if(particles_[k].p_.y_ == 0.0) particles_[k].p_.y_ = 1.0;
        if(particles_[k].p_.z_ == 0.0) particles_[k].p_.z_ = 1.0;

        if(dim_ < 3) particles_[k].p_.z_ = 0.0;
        if(dim_ < 2) particles_[k].p_.y_ = 0.0;
    }

    // Reset momentum of all atoms that are not mobile
    for(int i=0; i<N_; i++)
    {
        if(state_->atoms_[i]->r_.size() == 0) continue;
        if(!state_->atoms_[i]->isMobile_) particles_[i].p_.set(0.0, 0.0, 0.0);
    }

    // Perform velocity scaling to achive desired initial T
    double Tnow = calcTemp();
    double p_scale = sqrt(NH_T_.T_ / Tnow);
    for(int k=0; k<N_; k++)
    {
        particles_[k].p_.x_*= p_scale;
        particles_[k].p_.y_*= p_scale;
        particles_[k].p_.z_*= p_scale;
    }
}

void CSimulationBox::pbcAdd(double& pos, double L, double Lm)
{
    // For particle PBCs, add/subtract the
    // appropriate number of integer PBCs to
    // translate particle to inside of PBC
    pos-= double(int((pos + L)/Lm))*Lm;
}

void CSimulationBox::pbcWrap()
{
    C3DVector simBox = calcSimBox();
    double LmX = simBox.x_;
    double LmY = simBox.y_;
    double LmZ = simBox.z_;
    double LmaxX_2 = LmX / 2.0;
    double LmaxY_2 = LmY / 2.0;
    double LmaxZ_2 = LmZ / 2.0;

    for(int k=0; k<N_; k++)
    {
        if(particles_[k].r_.x_ >= LmaxX_2) pbcAdd(particles_[k].r_.x_,  LmaxX_2, LmX);
        if(particles_[k].r_.x_ < -LmaxX_2) pbcAdd(particles_[k].r_.x_, -LmaxX_2, LmX);

        if(dim_ < 2) continue;
        if(particles_[k].r_.y_ >= LmaxY_2) pbcAdd(particles_[k].r_.y_,  LmaxY_2, LmY);
        if(particles_[k].r_.y_ < -LmaxY_2) pbcAdd(particles_[k].r_.y_, -LmaxY_2, LmY);

        if(dim_ < 3) continue;
        if(particles_[k].r_.z_ >= LmaxZ_2) pbcAdd(particles_[k].r_.z_,  LmaxZ_2, LmZ);
        if(particles_[k].r_.z_ < -LmaxZ_2) pbcAdd(particles_[k].r_.z_, -LmaxZ_2, LmZ);
    }
}

double CSimulationBox::calcTemp()
{
    double T = 0.0;
    int N = (int)particles_.size();

    for(int k=0; k<N; k++)
    {
        C3DVector p = particles_[k].p_;
        T+= (p*p) / particles_[k].m_;
    }
    
    return (T / (dim_*double(N)));
}

double CSimulationBox::calcPress(const mthost_vector<CMDFFMatrices::CForces>& F) const
{
    double V = velVerlet_.getV(LmaxX_, LmaxY_, LmaxZ_, molDynConfig_.ensemble_);
    double sum = 0.0;
    
    for(int k=0; k<N_; k++)
    {
        double p2 = (particles_[k].p_*particles_[k].p_) / particles_[k].m_;
        double f = F[k].Fpi_*particles_[k].r_;
        sum+= (p2 + f);
    }
    
    return (1.0 / (3.0*V)) * sum;
}

std::string CSimulationBox::getAtomType(int index)
{
    return state_->atoms_[index]->getID();
}

double CSimulationBox::calcCurrentKineticEnergy()
{
    double K = 0.0;
    for(CParticle3D& particle : particles_)
    {
        K+= particle.p_*particle.p_ / particle.m_;
    }

    return 0.5*K;
}

END_CUDA_COMPATIBLE()
