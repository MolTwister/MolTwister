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
    void PBCWrap();
    double CalcTemp();
    double CalcPress(const mthost_vector<CMDFFMatrices::CForces>& F) const;
    double CalcV() { return VelVerlet.GetV(LmaxX_, LmaxY_, LmaxZ_, bNPTEnsemble); }
    double getLmaxX() { return LmaxX_; }
    double getLmaxY() { return LmaxY_; }
    double getLmaxZ() { return LmaxZ_; }
    bool isNPTEnsemble() { return bNPTEnsemble; }

    void NHTPropagator(CFct& Fct)
        { NH_T.Propagator(N, dim, dt, Fct); }
    void NHPPropagator(CFct& Fct)
        { if(bNPTEnsemble) NH_P.Propagator(N, dim, dt, Fct); }
    void VelVerPropagator(mthost_vector<CMDFFMatrices::CForces>& F)
        { VelVerlet.Propagator(N, dim, dt, LmaxX_, LmaxY_, LmaxZ_, aParticles, F, bNPTEnsemble); }
    mthost_vector<CMDFFMatrices::CForces> CalcParticleForces()
        { return VelVerlet.CalcParticleForces(dim, LmaxX_, LmaxY_, LmaxZ_, aParticles); }
    
private:
    void InitSystem();
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
    double dt;                                      // Time-step in reduced units
    EDim dim;                                       // System dimension
    bool bNegMomHalfWay;                            // If true, p=-p for all p, midway the simulation

private:
    double LmaxX_;                                   // Box side length in Å
    double LmaxY_;                                   // Box side length in Å
    double LmaxZ_;                                   // Box side length in Å
    bool bNPTEnsemble;                              // If true, NPT, else NVT
    CMolTwisterState* state_;
    SMolDynConfigStruct molDynConfig_;
};

END_CUDA_COMPATIBLE()
