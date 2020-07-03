#pragma once
#include "../../../Utilities/3DVector.h"
#include "../Integrators/Particle3D.h"
#include "../Integrators/NHChain.h"
#include "../Integrators/VelVerlet.h"
#include "../../../MolTwisterState.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CSimulationBox
{
public:
    enum EDim { dim1D = 1, dim2D = 2, dim3D = 3 };

public:
    CSimulationBox(CMolTwisterState* state, FILE* stdOut);
    
public:
    void InitSystem(int iM);
    void PBCWrap();
    double CalcTemp();
    double CalcPress(const mthost_vector<CMDFFMatrices::CForces>& F) const;
    double CalcV() { return VelVerlet.GetV(LmaxX, LmaxY, LmaxZ, bNPTEnsemble); }

    void NHTPropagator(CFct& Fct)
        { NH_T.Propagator(N, dim, dt, Fct); }
    void NHPPropagator(CFct& Fct)
        { if(bNPTEnsemble) NH_P.Propagator(N, dim, dt, Fct); }
    void VelVerPropagator(mthost_vector<CMDFFMatrices::CForces>& F)
        { VelVerlet.Propagator(N, dim, dt, LmaxX, LmaxY, LmaxZ, aParticles, F, bNPTEnsemble); }
    mthost_vector<CMDFFMatrices::CForces> CalcParticleForces()
        { return VelVerlet.CalcParticleForces(dim, LmaxX, LmaxY, LmaxZ, aParticles); }
    
private:
    void ResizeArrays();

    void copyPosFromState();
    void SetRandMom(int iFirst, int iTo);
    void PBCAdd(double& pos, double L, double Lm);

public:
    mthost_vector<CParticle3D> aParticles;          // Particle positions, velocities, masses etc.
    CNHChain NH_T;                                  // Nose-Hoover chain propagator temperature
    CNHChain NH_P;                                  // Nose-Hoover chain propagator pressure
    CVelVerlet VelVerlet;                           // Velocity verlet propagator
    int N;                                          // Number of particles
    double LmaxX;                                   // Box side length in Å
    double LmaxY;                                   // Box side length in Å
    double LmaxZ;                                   // Box side length in Å
    double dt;                                      // Time-step in reduced units
    EDim dim;                                       // System dimension
    bool bNegMomHalfWay;                            // If true, p=-p for all p, midway the simulation
    bool bNPTEnsemble;                              // If true, NPT, else NVT

private:
    CMolTwisterState* state_;
};

END_CUDA_COMPATIBLE()
