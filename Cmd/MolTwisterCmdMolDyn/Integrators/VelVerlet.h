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
#include <stdio.h>
#include <vector>
#include "Particle3D.h"
#include "Constants.h"
#include "../ForceFields/MDFFMatrices.h"
#include "../../../MolTwisterState.h"
#include "../../../Utilities/CUDAGeneralizations.h"
#include "../Config/MolDynConfigStruct.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFMatrices;
class CVelVerlet
{
public:
    CVelVerlet(CMolTwisterState* state, FILE* stdOut, double rCutoff, double dShell, double fCutoff);
    virtual ~CVelVerlet();

public:
    void propagator(int N, int dim, double dt, double LmaxX, double LmaxY, double LmaxZ, mthost_vector<CParticle3D>& particles, mthost_vector<CMDFFMatrices::CForces>& F, SMolDynConfigStruct::Ensemble ensemble, C3DVector& boxSizeOut);
    mthost_vector<CMDFFMatrices::CForces> calcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& particles);
    void setRandMom(double tau);
    void setNonBondedScaleFactors(float scale12, float scale13, float scale14, float scale1N);
    double getV(double LmaxX, double LmaxY, double LmaxZ, SMolDynConfigStruct::Ensemble ensemble) const;
    C3DVector getSimBox(double LmaxX, double LmaxY, double LmaxZ, SMolDynConfigStruct::Ensemble ensemble) const;
    void setVerboseOutput(bool verbose=true) { verboseOutput_ = verbose; }

private:
    void calcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& particles, mthost_vector<CMDFFMatrices::CForces>& F);
    double G_eps(int N, const mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>& F);
    void cutForces(mthost_vector<CMDFFMatrices::CForces>& F);
    void prop_p(int N, double dt, mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>& F);
    void prop_r(int N, double dt, mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>& F);

public:
    double Fcut_;
    double Rcut_;
    double P_;
    double V0_;
    double W_;
    double p_eps_;

private:
    double eps_;
    CMDFFMatrices* mdFFMatrices_;
    float scale12_;
    float scale13_;
    float scale14_;
    float scale1N_;
    bool verboseOutput_;
};

END_CUDA_COMPATIBLE()
