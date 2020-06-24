#pragma once
#include <stdio.h>
#include <vector>
#include "Particle3D.h"
#include "Constants.h"
#include "../ForceFields/MDForces.h"
#include "../../../MolTwisterState.h"

class CMDFFMatrices;
class CVelVerlet
{
public:
    CVelVerlet(CMolTwisterState* state, FILE* stdOut);
    virtual ~CVelVerlet();

public:
    void Propagator(int N, int dim, double dt, double Lmax, std::vector<CParticle3D>& aParticles, std::vector<CMDForces>& F, bool bNPT=false);
    std::vector<CMDForces> CalcParticleForces(int N, int dim, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles);
    void SetRandMom(double tau);
    double GetV(double Lmax, bool bNPT=false) const;

private:
    double G_eps(int N, const std::vector<CParticle3D>& aParticles, const std::vector<CMDForces>& F);
    void Prop_p(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<CMDForces>& F);
    void Prop_r(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<CMDForces>& F);

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
