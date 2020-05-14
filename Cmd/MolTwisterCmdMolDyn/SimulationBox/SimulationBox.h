#ifndef __ThesisMDTests__SimulationBox__
#define __ThesisMDTests__SimulationBox__

#include "../../../Utilities/3DVector.h"
#include "../Integrators/Particle3D.h"
#include "../ForceFields/ForceFields.h"
#include "../Integrators/NHChain.h"
#include "../Integrators/VelVerlet.h"

class CSimulationBox
{
public:
    enum ESystem { sys1DHarmBond = 1, sys1DExtHarmOsc = 2, sys3DExtHarmOsc = 3,
        sysLJCH4HighDens = 4, sysLJCH4NormDens = 5, sys1DFree = 6 };
    enum EDim { dim1D = 1, dim2D = 2, dim3D = 3 };

public:
    CSimulationBox();
    
public:
    void InitSystem(ESystem System, int iM);
    void PBCWrap();
    double CalcTemp();
    double CalcPress(vector<C3DVector>& aF);
    double CalcV() { return VelVerlet.GetV(Lmax, bNPTEnsemble); }

    void NHTPropagator(CFct& Fct)
        { NH_T.Propagator(N, dim, dt, Fct); }
    void NHPPropagator(CFct& Fct)
        { if(bNPTEnsemble) NH_P.Propagator(N, dim, dt, Fct); }
    void VelVerPropagator(vector<C3DVector>& aF, vector<C3DVector>& aFpi)
        { VelVerlet.Propagator(N, dim, dt, Lmax, aParticles, aF, aFpi, bNPTEnsemble); }
    C3DVector CalcParticleForce(int k, C3DVector& Fpi)
        { return VelVerlet.CalcParticleForce(k, N, dim, Lmax, Lmax, Lmax, aParticles, Fpi); }
    
private:
    void ResizeArrays();
    
public:
    vector<CParticle3D> aParticles;                 // Particle positions, velocities, masses etc.
    CNHChain NH_T;                                  // Nose-Hoover chain propagator temperature
    CNHChain NH_P;                                  // Nose-Hoover chain propagator pressure
    CVelVerlet VelVerlet;                           // Velocity verlet propagator
    int N;                                          // Number of particles
    double Lmax;                                    // Box side length in Å
    double dt;                                      // Time-step in reduced units
    EDim dim;                                       // System dimension
    bool bNegMomHalfWay;                            // If true, p=-p for all p, midway the simulation
    bool bNPTEnsemble;                              // If true, NPT, else NVT
    
private:
    void Init1DHarmBond();
    void Init1DExtHarmOsc();
    void Init3DExtHarmOsc();
    void InitLJCH4(double dBoxLen);
    void InitFreeParticles1D();
    void SetMasses(double m, int iFirst, int iTo);
    void SetRandPos(int iFirst, int iTo);
    void Set3DLatticePos(int iSize);
    void SetRandMom(int iFirst, int iTo);
    void PBCAdd(double& pos, double L, double Lm);
};

#endif
