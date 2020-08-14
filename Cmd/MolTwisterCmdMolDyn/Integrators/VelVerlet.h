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
    void Propagator(int N, int dim, double dt, double LmaxX, double LmaxY, double LmaxZ, mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDFFMatrices::CForces>& F, SMolDynConfigStruct::Ensemble ensemble);
    mthost_vector<CMDFFMatrices::CForces> CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles);
    void SetRandMom(double tau);
    double GetV(double LmaxX, double LmaxY, double LmaxZ, SMolDynConfigStruct::Ensemble ensemble) const;

private:
    void CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDFFMatrices::CForces>& F);
    double G_eps(int N, const mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>& F);
    void Prop_p(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>& F);
    void Prop_r(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>& F);

public:
    double Fcut;
    double P;
    double V0;
    double W;
    double p_eps;

private:
    double eps;
    CMDFFMatrices* mdFFMatrices_;
};

END_CUDA_COMPATIBLE()
