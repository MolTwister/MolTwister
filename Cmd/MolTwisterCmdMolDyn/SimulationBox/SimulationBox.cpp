#include "SimulationBox.h"
#include "../Integrators/Constants.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

CSimulationBox::CSimulationBox(CMolTwisterState* state, FILE* stdOut, SMolDynConfigStruct config)
    : VelVerlet(state, stdOut, config.cutoffRadius_, config.neighListShell_, config.cutoffForce_)
{
    state_ = state;
    molDynConfig_ = config;

    N = 100;
    C3DRect pbc = state->view3D_->getPBC();
    LmaxX_ = pbc.getWidthX(); // [Å]
    LmaxY_ = pbc.getWidthY(); // [Å]
    LmaxZ_ = pbc.getWidthZ(); // [Å]
    dt = 1.0 / Conv_t;   // [fs]
    dim = dim3D;
    bNegMomHalfWay = false;
    ResizeArrays();

    InitSystem();
}

void CSimulationBox::InitSystem()
{
    NH_T.M = molDynConfig_.temperatureNHChainLength_;
    NH_P.M = molDynConfig_.pressureNHChainLength_;

    NH_T.n = molDynConfig_.temperatureRESPASteps_;
    NH_P.n = molDynConfig_.pressureRESPASteps_;
    N = (int)state_->atoms_.size();
    NH_T.T = NH_P.T = molDynConfig_.temperature_ / Conv_T;
    dt = molDynConfig_.timeStep_ / Conv_t;
    NH_T.tau = molDynConfig_.temperatureRelaxationTime_ / Conv_t;
    NH_P.tau = molDynConfig_.pressureRelaxationTime_ / Conv_t;
    dim = dim3D;
    ResizeArrays();
    NH_T.SetRandNHPos();
    NH_P.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();
    NH_P.SetRandNHMom();

    VelVerlet.W = NH_P.T*NH_P.tau*NH_P.tau;
    VelVerlet.SetRandMom(NH_P.tau);
    VelVerlet.P = molDynConfig_.pressure_ / Conv_press;
    VelVerlet.V0 = LmaxX_*LmaxY_*LmaxZ_;

    copyPosFromState();
}

void CSimulationBox::ResizeArrays()
{
    NH_T.PrepareArrays(N, dim);
    NH_P.PrepareArrays(1, 1);
    aParticles.resize(N);
}

void CSimulationBox::copyPosFromState()
{
    for(int i=0; i<N; i++)
    {
        if(state_->atoms_[i]->r_.size() == 0) continue;
        aParticles[i].m = state_->atoms_[i]->m_;
        aParticles[i].x = state_->atoms_[i]->r_[0];
    }
}

void CSimulationBox::SetRandMom(int iFirst, int iTo)
{
    double          a = 2.0 / double(RAND_MAX);
    
    // Randomize p between 0 and 1 and prevent p=0
    for(int k=iFirst; k<iTo; k++)
    {
        aParticles[k].p.x_ = (a*double(rand()) - 1.0);
        aParticles[k].p.y_ = (a*double(rand()) - 1.0);
        aParticles[k].p.z_ = (a*double(rand()) - 1.0);
        
        if(aParticles[k].p.x_ == 0.0) aParticles[k].p.x_ = 1.0;
        if(aParticles[k].p.y_ == 0.0) aParticles[k].p.y_ = 1.0;
        if(aParticles[k].p.z_ == 0.0) aParticles[k].p.z_ = 1.0;

        if(dim < 3) aParticles[k].p.z_ = 0.0;
        if(dim < 2) aParticles[k].p.y_ = 0.0;
    }

    // Perform velocity scaling to achive desired initial T
    double Tnow = CalcTemp();
    double p_scale = sqrt(NH_T.T / Tnow);
    for(int k=iFirst; k<iTo; k++)
    {
        aParticles[k].p.x_*= p_scale;
        aParticles[k].p.y_*= p_scale;
        aParticles[k].p.z_*= p_scale;
    }
}

void CSimulationBox::PBCAdd(double& pos, double L, double Lm)
{
    // For particle PBCs, add/subtract the
    // appropriate number of integer PBCs to
    // translate particle to inside of PBC
    pos-= double(int((pos + L)/Lm))*Lm;
}

void CSimulationBox::PBCWrap()
{
    double Lm = pow(CalcV(), 1.0 / 3.0);
    double Lmax_2 = Lm / 2.0;

    for(int k=0; k<N; k++)
    {
        if(aParticles[k].x.x_ >= Lmax_2) PBCAdd(aParticles[k].x.x_,  Lmax_2, Lm);
        if(aParticles[k].x.x_ < -Lmax_2) PBCAdd(aParticles[k].x.x_, -Lmax_2, Lm);

        if(dim < 2) continue;
        if(aParticles[k].x.y_ >= Lmax_2) PBCAdd(aParticles[k].x.y_,  Lmax_2, Lm);
        if(aParticles[k].x.y_ < -Lmax_2) PBCAdd(aParticles[k].x.y_, -Lmax_2, Lm);

        if(dim < 3) continue;
        if(aParticles[k].x.z_ >= Lmax_2) PBCAdd(aParticles[k].x.z_,  Lmax_2, Lm);
        if(aParticles[k].x.z_ < -Lmax_2) PBCAdd(aParticles[k].x.z_, -Lmax_2, Lm);
    }
}

double CSimulationBox::CalcTemp()
{
    double  T = 0.0;
    int     N = (int)aParticles.size();

    for(int k=0; k<N; k++)
    {
        C3DVector p = aParticles[k].p;
        T+= (p*p) / aParticles[k].m;
    }
    
    return (T / (dim*double(N)));
}

double CSimulationBox::CalcPress(const mthost_vector<CMDFFMatrices::CForces>& F) const
{
    double  V = VelVerlet.GetV(LmaxX_, LmaxY_, LmaxZ_, molDynConfig_.ensemble_);
    double  sum = 0.0;
    
    for(int k=0; k<N; k++)
    {
        double p2 = (aParticles[k].p*aParticles[k].p) / aParticles[k].m;
        double f = F[k].Fpi_*aParticles[k].x;
        sum+= (p2 + f);
    }
    
    return (1.0 / (3.0*V)) * sum;
}

END_CUDA_COMPATIBLE()
