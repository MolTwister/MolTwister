#include "SimulationBox.h"
#include "../Integrators/Constants.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

CSimulationBox::CSimulationBox(CMolTwisterState* state, FILE* stdOut) : VelVerlet(state, stdOut)
{
    state_ = state;

    N = 100;
    LmaxX = 40.0;        // [Å]
    LmaxY = 40.0;        // [Å]
    LmaxZ = 40.0;        // [Å]
    dt = 1.0 / Conv_t;   // [fs]
    dim = dim3D;
    bNegMomHalfWay = false;
    bNPTEnsemble = false;
    ResizeArrays();
}
/*
void CSimulationBox::Init1DHarmBond()
{
    ////////////////////////////////////////
    // Two CH4 molecules in harmonic bond
    // K = 0.6 kJ/(mol*Å^2), r0 = 1.0Å
    ////////////////////////////////////////

    NH_T.n = 4;
    N = 2;
    Lmax = 40.0;
    NH_T.T = 298.0 / Conv_T;
    dt = 1.0 / Conv_t;
    NH_T.tau = 20.0*dt;
    dim = dim1D;
    ResizeArrays();
    SetMasses(16.04, 0, N);
    NH_T.SetRandNHPos();
    NH_T.SetRandNHMom();

    aParticles[0].p.x_ = 0.0;
    aParticles[1].p.x_ = 0.0;
    aParticles[0].x.x_ =  1.0;
    aParticles[1].x.x_ = -1.0;
    
    VelVerlet.aFHarmBond.resize(1);
    VelVerlet.aFHarmBond[0].Set(100.0, 1.0, 0, 1);
}

void CSimulationBox::Init1DExtHarmOsc()
{
    ////////////////////////////////////////
    // Free CH4 gas trapped in harm. osc.
    // potential at center of sim. line
    // K = 0.6 kJ/(mol*Å^2)
    ////////////////////////////////////////
    
    NH_T.n = 4;
    N = 1;
    Lmax = 40.0;
    NH_T.T = 298.0 / Conv_T;
    dt = 1.0 / Conv_t;
    NH_T.tau = 20.0*dt;
    dim = dim1D;
    ResizeArrays();
    SetMasses(16.04, 0, N);
    SetRandPos(0, N);
    NH_T.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();

    C3DVector r0 = C3DVector(0.0, 0.0, 0.0);
    for(int k=0; k<N; k++) VelVerlet.aFExternal[k].Harm.Set(0.6, r0);
}

void CSimulationBox::Init3DExtHarmOsc()
{
    ////////////////////////////////////////
    // Free CH4 gas trapped in harm. osc.
    // potential at center of sim. box
    // K = 0.6 kJ/(mol*Å^2)
    ////////////////////////////////////////

    NH_T.n = 4;
    N = 100;
    Lmax = 40.0;
    NH_T.T = 298.0 / Conv_T;
    dt = 1.0 / Conv_t;
    NH_T.tau = 20.0*dt;
    dim = dim3D;
    ResizeArrays();
    SetMasses(16.04, 0, N);
    SetRandPos(0, N);
    NH_T.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();
    
    C3DVector r0 = C3DVector(0.0, 0.0, 0.0);
    for(int k=0; k<N; k++) VelVerlet.aFExternal[k].Harm.Set(0.6, r0);
}

void CSimulationBox::InitLJCH4(double dBoxLen)
{
    ////////////////////////////////////////
    // OPLS-UA LJ CH4 gas in cubic system
    // sigma=3.730 Å, epsilon=1.2301 kJ/mol
    ////////////////////////////////////////

    int iLattSize = 5;

    NH_T.n = NH_P.n = 4;
    N = iLattSize*iLattSize*iLattSize;
    NH_T.T = NH_P.T = 298.0 / Conv_T;
    dt = 0.1 / Conv_t;
    double tauP = 5000.0*dt;
    NH_T.tau = 20.0*dt;
    NH_P.tau = tauP;
    Lmax = dBoxLen;
    dim = dim3D;
    ResizeArrays();
    SetMasses(16.04, 0, N);
    Set3DLatticePos(iLattSize);
    NH_T.SetRandNHPos();
    NH_P.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();
    NH_P.SetRandNHMom();

    VelVerlet.W = NH_P.T*tauP*tauP;
    VelVerlet.SetRandMom(tauP);
    VelVerlet.P = 10.0 / Conv_press;
    VelVerlet.V0 = Lmax*Lmax*Lmax;
    
    for(int k1=0; k1<N; k1++)
    {
        for(int k2=0; k2<N; k2++)
        {
            if(k1 != k2)
            {
                VelVerlet.aFNonBonded[k1][k2].LJ.Set(1.2301, 3.730);
            }
        }
    }
}

void CSimulationBox::InitFreeParticles1D()
{
    ////////////////////////////////////////
    // Syste of free CH4 gas particles in
    // cubic box with sides of 40Å
    ////////////////////////////////////////

    NH_T.n = 4;
    N = 1;
    Lmax = 40.0;
    NH_T.T = 298.0 / Conv_T;
    dt = 1.0 / Conv_t;
    NH_T.tau = 20.0*dt;
    dim = dim1D;
    bNegMomHalfWay = true;
    ResizeArrays();
    SetMasses(16.04, 0, N);
    SetRandPos(0, N);
    NH_T.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();
}
*/
void CSimulationBox::InitSystem(int iM)
{
    NH_T.M = iM;
    NH_P.M = iM;
    /*
    if(System == sys1DHarmBond) Init1DHarmBond();
    else if(System == sys1DExtHarmOsc) Init1DExtHarmOsc();
    else if(System == sys3DExtHarmOsc) Init3DExtHarmOsc();
    else if(System == sysLJCH4HighDens) InitLJCH4(35.0); // 77.7kg/m^3
    else if(System == sysLJCH4NormDens) InitLJCH4(172.0); // 0.65kg/m^3
    else if(System == sys1DFree) InitFreeParticles1D();*/

    ////////////////////////////////////////
    // OPLS-UA LJ CH4 gas in cubic system
    // sigma=3.730 Å, epsilon=1.2301 kJ/mol
    ////////////////////////////////////////

    // Configure MD parameters
    // :TODO: Need to set
    //          * correct box size from state (and not just cubic)
    //          * temperature from settings
    //          * Pressure from settings
    //          * NH-parameters from settings
    const double dBoxLen = 35.0;

    NH_T.n = NH_P.n = 4;
    N = (int)state_->atoms_.size();
    NH_T.T = NH_P.T = 298.0 / Conv_T;
    dt = 0.1 / Conv_t;
    double tauP = 5000.0*dt;
    NH_T.tau = 20.0*dt;
    NH_P.tau = tauP;
    LmaxX = dBoxLen;
    LmaxY = dBoxLen;
    LmaxZ = dBoxLen;
    dim = dim3D;
    ResizeArrays();
    NH_T.SetRandNHPos();
    NH_P.SetRandNHPos();
    SetRandMom(0, N);
    NH_T.SetRandNHMom();
    NH_P.SetRandNHMom();

    VelVerlet.W = NH_P.T*tauP*tauP;
    VelVerlet.SetRandMom(tauP);
    VelVerlet.P = 10.0 / Conv_press;
    VelVerlet.V0 = LmaxX*LmaxY*LmaxZ;

    copyPosFromState();
}

void CSimulationBox::ResizeArrays()
{
    NH_T.PrepareArrays(N, dim);
    NH_P.PrepareArrays(1, 1);
    aParticles.resize(N);
//    VelVerlet.aFExternal.resize(N);
//    VelVerlet.aFNonBonded.resize(N);
//    for(int k=0; k<N; k++) VelVerlet.aFNonBonded[k].resize(N);
}
/*
void CSimulationBox::SetMasses(double m, int iFirst, int iTo)
{
    for(int k=iFirst; k<iTo; k++) aParticles[k].m = m;
}
*/
/*
void CSimulationBox::SetRandPos(int iFirst, int iTo)
{
    double a = 2.0 / double(RAND_MAX);
    double b = Lmax / 2.0;
    
    for(int k=iFirst; k<iTo; k++)
    {
        aParticles[k].x.x_ = (a * double(rand()) - 1.0) * b;
        aParticles[k].x.y_ = (a * double(rand()) - 1.0) * b;
        aParticles[k].x.z_ = (a * double(rand()) - 1.0) * b;
        
        if(dim < 3) aParticles[k].x.z_ = 0.0;
        if(dim < 2) aParticles[k].x.y_ = 0.0;
    }
}
*/
void CSimulationBox::copyPosFromState()
{
    for(int i=0; i<N; i++)
    {
        if(state_->atoms_[i]->r_.size() == 0) continue;
        aParticles[i].m = state_->atoms_[i]->m_;
        aParticles[i].x = state_->atoms_[i]->r_[0];
    }
}

/*
void CSimulationBox::Set3DLatticePos(int iSize)
{
    int     iIndex = 0;
    double  a = double(iSize);
    double  b = 0.5 * (1.0 - a) / a;
    
    for(int ix=0; ix<iSize; ix++)
    {
        for(int iy=0; iy<iSize; iy++)
        {
            for(int iz=0; iz<iSize; iz++)
            {
                aParticles[iIndex].x.x_ = Lmax * (double(ix) / a + b);
                aParticles[iIndex].x.y_ = Lmax * (double(iy) / a + b);
                aParticles[iIndex].x.z_ = Lmax * (double(iz) / a + b);
                
                iIndex++;
            }
        }
    }
}
*/
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
    double  V = VelVerlet.GetV(LmaxX, LmaxY, LmaxZ, bNPTEnsemble);
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
