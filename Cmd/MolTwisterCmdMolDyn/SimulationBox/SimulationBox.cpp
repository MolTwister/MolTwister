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

    N_ = 100;
    C3DRect pbc = state->view3D_->getPBC();
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
    setRandMom(0, N_);
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
        particles_[i].x_ = state_->atoms_[i]->r_[0];
    }
}

void CSimulationBox::setRandMom(int first, int to)
{
    double  a = 2.0 / double(RAND_MAX);
    
    // Randomize p between 0 and 1 and prevent p=0
    for(int k=first; k<to; k++)
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

    // Perform velocity scaling to achive desired initial T
    double Tnow = calcTemp();
    double p_scale = sqrt(NH_T_.T_ / Tnow);
    for(int k=first; k<to; k++)
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
    double Lm = pow(calcV(), 1.0 / 3.0);
    double Lmax_2 = Lm / 2.0;

    for(int k=0; k<N_; k++)
    {
        if(particles_[k].x_.x_ >= Lmax_2) pbcAdd(particles_[k].x_.x_,  Lmax_2, Lm);
        if(particles_[k].x_.x_ < -Lmax_2) pbcAdd(particles_[k].x_.x_, -Lmax_2, Lm);

        if(dim_ < 2) continue;
        if(particles_[k].x_.y_ >= Lmax_2) pbcAdd(particles_[k].x_.y_,  Lmax_2, Lm);
        if(particles_[k].x_.y_ < -Lmax_2) pbcAdd(particles_[k].x_.y_, -Lmax_2, Lm);

        if(dim_ < 3) continue;
        if(particles_[k].x_.z_ >= Lmax_2) pbcAdd(particles_[k].x_.z_,  Lmax_2, Lm);
        if(particles_[k].x_.z_ < -Lmax_2) pbcAdd(particles_[k].x_.z_, -Lmax_2, Lm);
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
        double f = F[k].Fpi_*particles_[k].x_;
        sum+= (p2 + f);
    }
    
    return (1.0 / (3.0*V)) * sum;
}

std::string CSimulationBox::getAtomType(int index)
{
    return state_->atoms_[index]->getID();
}

END_CUDA_COMPATIBLE()
