#include "VelVerlet.h"
#include "Math.h"
#include <float.h>
#include <math.h>
#include <algorithm>
#include "../ForceFields/FunctorCalcForce.h"
#include "../MDLoop/Printf.h"

BEGIN_CUDA_COMPATIBLE()

CVelVerlet::CVelVerlet(CMolTwisterState* state, FILE* stdOut, double rCutoff, double dShell, double fCutoff)
{
    verboseOutput_ = false;

    P_ = 1.0 / Conv_press;
    p_eps_ = 1.0;
    eps_ = 0.0;
    W_ = 1.0;
    Fcut_ = fCutoff;
    Rcut_ = rCutoff;

    scale12_ = 0.0f;
    scale13_ = 0.0f;
    scale14_ = 0.5f;
    scale1N_ = 0.0f;

    mdFFMatrices_ = new CMDFFMatrices(state, stdOut, (float)rCutoff, (float)dShell);
}

CVelVerlet::~CVelVerlet()
{
    if(mdFFMatrices_) delete mdFFMatrices_;
}

void CVelVerlet::propagator(int N, int dim, double dt, double LmaxX, double LmaxY, double LmaxZ, mthost_vector<CParticle3D>& particles, mthost_vector<CMDFFMatrices::CForces>& F, SMolDynConfigStruct::Ensemble ensemble, C3DVector& boxSizeOut)
{
    // This is the implementation for NPT
    if(ensemble == SMolDynConfigStruct::ensembleNPT)
    {
        p_eps_+= ((dt / 2.0) * G_eps(N, particles, F)); // Step 3.1
        prop_p(N, dt, particles, F);                    // Step 3.2
        prop_r(N, dt, particles, F);                    // Step 3.3
        eps_+= (dt * p_eps_ / W_);                      // Step 3.4

        boxSizeOut = getSimBox(LmaxX, LmaxY, LmaxZ, ensemble);

        calcParticleForces(dim, boxSizeOut.x_, boxSizeOut.y_, boxSizeOut.z_, particles, F);

        prop_p(N, dt, particles, F);                     // Step 3.5
        p_eps_+= ((dt / 2.0) * G_eps(N, particles, F));  // Step 3.6
    }

    // This is the same implementation for both NVT and NVE
    else
    {
        double dt_2 = dt / 2.0;
        for(int k=0; k<N; k++)
        {
            double m_k = particles[k].m_;
            particles[k].p_+= F[k].F_*dt_2;
            particles[k].r_+= particles[k].p_*(dt / m_k);
        }

        boxSizeOut = getSimBox(LmaxX, LmaxY, LmaxZ, ensemble);

        calcParticleForces(dim, LmaxX, LmaxY, LmaxZ, particles, F);
        if(Fcut_ > 0) cutForces(F);

        for(int k=0; k<N; k++)
        {
            particles[k].p_+= F[k].F_*dt_2;
        }
    }
}

void CVelVerlet::cutForces(mthost_vector<CMDFFMatrices::CForces>& F)
{
    for(size_t i=0; i<F.size(); i++)
    {
        if(fabs(F[i].F_.x_) > Fcut_)
        {
            double sign = ((F[i].F_.x_ >= 0.0) ? 1.0 : -1.0);
            double cutF = sign * Fcut_;
            if(verboseOutput_) COut::printf("Warning: x-forces of atom %i will be cut from %g to %g\r\n", (int)i, F[i].F_.x_, cutF);
            F[i].Fpi_.x_ = F[i].F_.x_ = cutF;
        }
        if(fabs(F[i].F_.y_) > Fcut_)
        {
            double sign = ((F[i].F_.y_ >= 0.0) ? 1.0 : -1.0);
            double cutF = sign * Fcut_;
            if(verboseOutput_) COut::printf("Warning: y-forces of atom %i will be cut from %g to %g\r\n", (int)i, F[i].F_.y_, cutF);
            F[i].Fpi_.y_ = F[i].F_.y_ = cutF;
        }
        if(fabs(F[i].F_.z_) > Fcut_)
        {
            double sign = ((F[i].F_.z_ >= 0.0) ? 1.0 : -1.0);
            double cutF = sign * Fcut_;
            if(verboseOutput_) COut::printf("Warning: z-forces of atom %i will be cut from %g to %g\r\n", (int)i, F[i].F_.z_, cutF);
            F[i].Fpi_.z_ = F[i].F_.z_ = cutF;
        }
    }
}

void CVelVerlet::prop_p(int N, double dt, mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>& F)
{
    double u_eps = p_eps_ / W_;
    double alpha = (1.0 + 1.0/double(N));
    double parm = alpha*u_eps*dt / 4.0;
    double coeff = CMt::exp(-parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::sinhXoverX(parm); // sinh(parm) / parm;
    for(int k=0; k<N; k++)
    {
        particles[k].p_ = particles[k].p_*coeff2 + F[k].F_*((dt/2.0)*coeff*coeff3);
    }
}

void CVelVerlet::prop_r(int N, double dt, mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>&)
{
    double u_eps = p_eps_ / W_;
    double parm = u_eps*dt / 2.0;
    double coeff = CMt::exp(parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::sinhXoverX(parm) * coeff;
    for(int k=0; k<N; k++)
    {
        C3DVector u_k = particles[k].p_ * (1.0 / particles[k].m_);
        particles[k].r_ = particles[k].r_*coeff2 + u_k*(dt*coeff3);
    }
}

mthost_vector<CMDFFMatrices::CForces> CVelVerlet::calcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& particles)
{
    mthost_vector<CMDFFMatrices::CForces> F(particles.size());
    calcParticleForces(dim, Lx, Ly, Lz, particles, F);

    return F;
}

double CVelVerlet::G_eps(int N, const mthost_vector<CParticle3D>& particles, const mthost_vector<CMDFFMatrices::CForces>& F)
{
    double V = V0_ * exp(3.0 * eps_);

    double sum_p = 0.0;
    double sum_f = 0.0;
    for(int k=0; k<N; k++)
    {
        double p2 = particles[k].p_ * particles[k].p_;
        sum_p+= (p2 / particles[k].m_);
        sum_f+= (F[k].Fpi_ * particles[k].r_);
    }

    return (1.0 + 1.0 / double(N))*sum_p + sum_f - 3.0*P_*V;
}

void CVelVerlet::setRandMom(double tau)
{
    double a = 2.0 / double(RAND_MAX);

    p_eps_ = (a*double(rand()) - 1.0) * (W_ / tau);
    if(p_eps_ == 0.0) p_eps_ = (W_ / tau);
}

void CVelVerlet::setNonBondedScaleFactors(float scale12, float scale13, float scale14, float scale1N)
{
    scale12_ = scale12;
    scale13_ = scale13;
    scale14_ = scale14;
    scale1N_ = scale1N;
}

double CVelVerlet::getV(double LmaxX, double LmaxY, double LmaxZ, SMolDynConfigStruct::Ensemble ensemble) const
{
    if(ensemble == SMolDynConfigStruct::ensembleNPT) return V0_ * exp(3.0 * eps_);

    return LmaxX*LmaxY*LmaxZ;
}

C3DVector CVelVerlet::getSimBox(double LmaxX, double LmaxY, double LmaxZ, SMolDynConfigStruct::Ensemble ensemble) const
{
    if(ensemble == SMolDynConfigStruct::ensembleNPT)
    {
        double Vmax = LmaxX * LmaxY * LmaxZ;
        double etaCube = getV(LmaxX, LmaxY, LmaxZ, SMolDynConfigStruct::ensembleNPT) / Vmax;
        etaCube = pow(etaCube, 1.0/3.0);

        return C3DVector(etaCube * LmaxX, etaCube * LmaxY, etaCube * LmaxZ);
    }

    return C3DVector(LmaxX, LmaxY, LmaxZ);
}

void CVelVerlet::calcParticleForces(int dim, double Lx, double Ly, double Lz, const mthost_vector<CParticle3D>& particles, mthost_vector<CMDFFMatrices::CForces>& F)
{
    mdFFMatrices_->updateAtomList(particles);
    mdFFMatrices_->genNeighList((float)Lx, (float)Ly, (float)Lz);
    CFunctorCalcForce calcForce(dim, (float)Lx, (float)Ly, (float)Lz, (float)Fcut_, (float)Rcut_, scale12_, scale13_, scale14_, scale1N_);
    calcForce.setForceFieldMatrices(*mdFFMatrices_);
    mttransform(EXEC_POLICY mdFFMatrices_->devAtomList_.begin(), mdFFMatrices_->devAtomList_.end(), mdFFMatrices_->devForcesList_.begin(), calcForce);
    mtcudaDeviceSynchronize();
    F = mdFFMatrices_->devForcesList_;
}

END_CUDA_COMPATIBLE()
