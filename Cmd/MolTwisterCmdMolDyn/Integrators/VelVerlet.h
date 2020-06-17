#pragma once
#include <stdio.h>
#include <vector>
#include "Particle3D.h"
#include "Constants.h"
#include "../../../Utilities/3DVector.h"
#include "../ForceFields/ForceFields.h"
#include "../../../MolTwisterState.h"

class CDevForces
{
public:
    CDevForces() {}
    CDevForces(const C3DVector& F, const C3DVector& Fpi) { F_ = F; Fpi_ = Fpi; }

public:
    C3DVector F_;
    C3DVector Fpi_;
};

class CDevForceFieldMatrices;
class CVelVerlet
{
public:
    CVelVerlet(CMolTwisterState* state, FILE* stdOut);
    virtual ~CVelVerlet();

public:
    void Propagator(int N, int dim, double dt, double Lmax, std::vector<CParticle3D>& aParticles, std::vector<CDevForces>& F, bool bNPT=false);
    std::vector<CDevForces> CalcParticleForces(int N, int dim, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles);
    void SetRandMom(double tau);
    double GetV(double Lmax, bool bNPT=false) const;

private:
    double G_eps(int N, const std::vector<CParticle3D>& aParticles, const std::vector<CDevForces>& F);
    void Prop_p(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<CDevForces>& F);
    void Prop_r(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<CDevForces>& F);
//    void StoreMaxF(const C3DVector& F);
//    void PrintDebugInfoAtCutForces(int k, int N, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles);

public:
    double Fcut;
    double P;
    double V0;
    double W;
    double p_eps;
    double eps;
    double tau;
//    std::vector<std::vector<CForceNonBonded>> aFNonBonded;  // Matrix of non-bonded forces between particles
//    std::vector<CForceExternal> aFExternal;                 // External forces assigned to each particle
//    std::vector<CForceHarmBond> aFHarmBond;                 // Harmonic bonds between particles

private:
    CDevForceFieldMatrices* devVelVerlet_;
};
