#pragma once
#include <stdio.h>
#include <vector>
#include "Particle3D.h"
#include "Constants.h"
#include "../../../Utilities/3DVector.h"
#include "../ForceFields/ForceFields.h"

class CVelVerlet
{
public:
    CVelVerlet() { P = 1.0 / Conv_press; p_eps = 1.0; eps = 0.0; W = 1.0; m_dLastFMax = 0.0; Fcut = 1000.0; m_bCutF = false; tau = 1.0; }
    
public:
    void Propagator(int N, int dim, double dt, double Lmax, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF, std::vector<C3DVector> &aFpi, bool bNPT=false);
    C3DVector CalcParticleForce(int k, int N, int dim, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles, C3DVector& Fpi);
    void SetRandMom(double tau);
    double GetV(double Lmax, bool bNPT=false);
    double GetMaxF() { return m_dLastFMax; }
    void PrintCutMsgAndReset();
    
private:
    double G_eps(int N, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF);
    void Prop_p(int N, double dt, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF);
    void Prop_r(int N, double dt, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF);
    void StoreMaxF(C3DVector& F);
    void PrintDebugInfoAtCutForces(int k, int N, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles);

public:
    double Fcut;
    double P;
    double V0;
    double W;
    double p_eps;
    double eps;
    double tau;
    std::vector<std::vector<CForceNonBonded>> aFNonBonded;  // Matrix of non-bonded forces between particles
    std::vector<CForceExternal> aFExternal;                 // External forces assigned to each particle
    std::vector<CForceHarmBond> aFHarmBond;                 // Harmonic bonds between particles
    
private:
    double m_dLastFMax;
    bool m_bCutF;
};
