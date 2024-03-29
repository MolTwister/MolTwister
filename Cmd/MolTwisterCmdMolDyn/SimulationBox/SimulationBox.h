//
// Copyright (C) 2023 Richard Olsen.
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

#pragma once
#include "../../../Utilities/3DVector.h"
#include "../Integrators/Particle3D.h"
#include "../Integrators/NHChain.h"
#include "../Integrators/VelVerlet.h"
#include "../../../MolTwisterState.h"
#include "../../../Utilities/CUDAGeneralizations.h"
#include "../Config/MolDynConfigStruct.h"

BEGIN_CUDA_COMPATIBLE()

class CSimulationBox
{
public:
    enum EDim { dim1D = 1, dim2D = 2, dim3D = 3 };

public:
    CSimulationBox(CMolTwisterState* state, FILE* stdOut, SMolDynConfigStruct config);
    
public:
    void pbcWrap();
    double calcTemp();
    double calcPress(const mthost_vector<CMDFFMatrices::CForces>& F) const;
    double calcV() { return velVerlet_.getV(LmaxX_, LmaxY_, LmaxZ_, molDynConfig_.ensemble_); }
    C3DVector calcSimBox() { return velVerlet_.getSimBox(LmaxX_, LmaxY_, LmaxZ_, molDynConfig_.ensemble_); }
    double getLmaxX() { return LmaxX_; }
    double getLmaxY() { return LmaxY_; }
    double getLmaxZ() { return LmaxZ_; }
    SMolDynConfigStruct::Ensemble getEnsemble() { return molDynConfig_.ensemble_; }
    double getTimeStep() { return molDynConfig_.timeStep_; }
    int getOutputStride() { return molDynConfig_.outputStride_; }
    int getNumTimeSteps() { return molDynConfig_.numberOfTimeSteps_; }
    double getGradDescentLearningRate() { return molDynConfig_.gradientDescentLearningRate_; }

    void NHTPropagator(CFct& fct)
        { if(molDynConfig_.ensemble_ == SMolDynConfigStruct::ensembleNPT || molDynConfig_.ensemble_ == SMolDynConfigStruct::ensembleNVT) NH_T_.propagator(N_, dim_, dt_, fct); }
    void NHPPropagator(CFct& fct)
        { if(molDynConfig_.ensemble_ == SMolDynConfigStruct::ensembleNPT) NH_P_.propagator(N_, dim_, dt_, fct); }
    void velVerPropagator(mthost_vector<CMDFFMatrices::CForces>& F, C3DVector& boxSizeOut)
        { velVerlet_.propagator(N_, dim_, dt_, LmaxX_, LmaxY_, LmaxZ_, particles_, F, molDynConfig_.ensemble_, boxSizeOut); }
    void gradientDescentStep(mthost_vector<CMDFFMatrices::CForces>& F);
    mthost_vector<CMDFFMatrices::CForces> calcParticleForces()
        { return velVerlet_.calcParticleForces(dim_, LmaxX_, LmaxY_, LmaxZ_, particles_); }
    std::string getAtomType(int index);
    double calcCurrentKineticEnergy();
    double calcUtot(const mthost_vector<CMDFFMatrices::CForces>& F) const;

private:
    void initSystem();
    void resizeArrays();

    void copyPosFromState();
    void setRandMom();
    void pbcAdd(double& pos, double L, double Lm);

public:
    mthost_vector<CParticle3D> particles_;           // Particle positions, velocities, masses etc.
    CNHChain NH_T_;                                  // Nose-Hoover chain propagator temperature
    CNHChain NH_P_;                                  // Nose-Hoover chain propagator pressure
    CVelVerlet velVerlet_;                           // Velocity verlet propagator
    int N_;                                          // Number of particles
    double dt_;                                      // Time-step in reduced units
    EDim dim_;                                       // System dimension

private:
    double LmaxX_;                                   // Box side length in AA
    double LmaxY_;                                   // Box side length in AA
    double LmaxZ_;                                   // Box side length in AA
    CMolTwisterState* state_;
    SMolDynConfigStruct molDynConfig_;
};

END_CUDA_COMPATIBLE()
