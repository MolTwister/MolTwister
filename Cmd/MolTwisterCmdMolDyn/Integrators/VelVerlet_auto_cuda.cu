#include "VelVerlet.h"
#include "Math.h"
#include <float.h>
#include <math.h>
#include <algorithm>
#include "../ForceFields/FunctorCalcForce.h"

BEGIN_CUDA_COMPATIBLE()

CVelVerlet::CVelVerlet(CMolTwisterState* state, FILE* stdOut)
{
    P = 1.0 / Conv_press;
    p_eps = 1.0;
    eps = 0.0;
    W = 1.0;
    Fcut = 1000.0;
    tau = 1.0;

    const float rCutoff = 10.0f;
    const float dShell = 2.0f;
    mdFFMatrices_ = new CMDFFMatrices(state, stdOut, rCutoff, dShell);
}

CVelVerlet::~CVelVerlet()
{
    if(mdFFMatrices_) delete mdFFMatrices_;
}

void CVelVerlet::Propagator(int N, int dim, double dt, double Lmax, mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDFFMatrices::CForces>& F, bool bNPT)
{
    if(!bNPT)
    {
        double dt_2 = dt / 2.0;
        for(int k=0; k<N; k++)
        {
            double m_k = aParticles[k].m;
            aParticles[k].p+= F[k].F_*dt_2;
            aParticles[k].x+= aParticles[k].p*(dt / m_k);
        }

        CalcParticleForces(dim, Lmax, Lmax, Lmax, aParticles, F);

        for(int k=0; k<N; k++)
        {
            aParticles[k].p+= F[k].F_*dt_2;
        }
    }

    else
    {
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, F));  // Step 3.1
        Prop_p(N, dt, aParticles, F);                    // Step 3.2
        Prop_r(N, dt, aParticles, F);                    // Step 3.3
        eps+= (dt * p_eps / W);                          // Step 3.4
        double Lm = pow(GetV(Lmax, true), 1.0/3.0);

        CalcParticleForces(dim, Lm, Lm, Lm, aParticles, F);

        Prop_p(N, dt, aParticles, F);                    // Step 3.5
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, F));  // Step 3.6
    }
}

void CVelVerlet::Prop_p(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>& F)
{
    double u_eps = p_eps / W;
    double alpha = (1.0 + 1.0/double(N));
    double parm = alpha*u_eps*dt / 4.0;
    double coeff = CMt::Exp(-parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm); // sinh(parm) / parm;
    for(int k=0; k<N; k++)
    {
        aParticles[k].p = aParticles[k].p*coeff2 + F[k].F_*((dt/2.0)*coeff*coeff3);
    }
}

void CVelVerlet::Prop_r(int N, double dt, mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>&)
{
    double u_eps = p_eps / W;
    double parm = u_eps*dt / 2.0;
    double coeff = CMt::Exp(parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm) * coeff;
    for(int k=0; k<N; k++)
    {
        C3DVector u_k = aParticles[k].p * (1.0 / aParticles[k].m);
        aParticles[k].x = aParticles[k].x*coeff2 + u_k*(dt*coeff3);
    }
}

mthost_vector<CMDFFMatrices::CForces> CVelVerlet::CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles)
{
    mthost_vector<CMDFFMatrices::CForces> F(aParticles.size());
    CalcParticleForces(dim, Lx, Ly, Lz, aParticles, F);

    return F;
}

double CVelVerlet::G_eps(int N, const mthost_vector<CParticle3D>& aParticles, const mthost_vector<CMDFFMatrices::CForces>& F)
{
    double V = V0 * exp(3.0 * eps);

    double sum_p = 0.0;
    double sum_f = 0.0;
    for(int k=0; k<N; k++)
    {
        double p2 = aParticles[k].p * aParticles[k].p;
        sum_p+= (p2 / aParticles[k].m);
        sum_f+= (F[k].Fpi_ * aParticles[k].x);
    }

    return (1.0 + 1.0 / double(N))*sum_p + sum_f - 3.0*P*V;
}

void CVelVerlet::SetRandMom(double tau)
{
    double a = 2.0 / double(RAND_MAX);

    p_eps = (a*double(rand()) - 1.0) * (W / tau);
    if(p_eps == 0.0) p_eps = (W / tau);
}

double CVelVerlet::GetV(double Lmax, bool bNPT) const
{
    if(!bNPT) return Lmax*Lmax*Lmax;

    return V0 * exp(3.0 * eps);
}

void CVelVerlet::CalcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& aParticles, mthost_vector<CMDFFMatrices::CForces>& F)
{
    mdFFMatrices_->updateAtomList(aParticles);
    CFunctorCalcForce calcForce(dim, Lx, Ly, Lz, Fcut);
    calcForce.setForceFieldMatrices(*mdFFMatrices_);
    mttransform(EXEC_POLICY mdFFMatrices_->devAtomList_.begin(), mdFFMatrices_->devAtomList_.end(), mdFFMatrices_->devForcesList_.begin(), calcForce);
    F = mdFFMatrices_->devForcesList_;
}

END_CUDA_COMPATIBLE()
