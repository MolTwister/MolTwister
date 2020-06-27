#pragma once
#include <stdio.h>
#include <vector>
#include "Particle3D.h"
#include "Constants.h"
#include "../ForceFields/MDForces.h"
#include "../../../MolTwisterState.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFMatrices;
class CVelVerlet
{
public:
    CVelVerlet(CMolTwisterState* state, FILE* stdOut);
    virtual ~CVelVerlet();

public:
    void Propagator(int N, int dim, double dt, double Lmax, mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDForces>& F, bool bNPT=false);
    mthost_vector<CMDForces> CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles);
    void SetRandMom(double tau);
    double GetV(double Lmax, bool bNPT=false) const;

private:
    void CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDForces>& F);
    double G_eps(int N, const mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDForces>& F);
    void Prop_p(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDForces>& F);
    void Prop_r(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDForces>& F);

public:
    double Fcut;
    double P;
    double V0;
    double W;
    double p_eps;
    double eps;
    double tau;

private:
    CMDFFMatrices* mdFFMatrices_;
};

END_CUDA_COMPATIBLE()
