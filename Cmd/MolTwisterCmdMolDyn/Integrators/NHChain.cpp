#include "NHChain.h"
#include "Constants.h"
#include "Math.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

BEGIN_CUDA_COMPATIBLE()

CNHChain::CNHChain()
{
    n_sy_ = 3;
    n_ = 4;
    M_ = 8;
    T_ = 298.0 / Conv_T;  // [K]
    tau_ = 20.0 / Conv_t;

    getSuzukiYoshida(w_);
}

void CNHChain::getSuzukiYoshida(double* w)
{
    const double w0 = 1.0 / (1.0 - pow(2.0, 1.0/3.0));
    
    w[0] = w0;
    w[1] = 1.0 - 2.0*w0;
    w[2] = w0;
}

void CNHChain::propagator(int N, int dim, double dt, CFct& f)
{
    double  beta = 1.0 / T_;
    double  coeff;
    double  pEtaCutoff = 15.0;
    
    
    if((M_ > 0) && (p_eta_[0] > pEtaCutoff))
    {
        printf("Warning! cutoff applied to NH chain (p_eta[0]=%g -> %g)...\r\n", p_eta_[0], pEtaCutoff);
        p_eta_[0] = 10.0;
    }
    
    for(int mu=0; mu<n_sy_; mu++)
    {
        double Dt = w_[mu] * dt / double(n_);
        for(int i=0; i<n_; i++)
        {
            // Step 1.1
            p_eta_[M_-1]+= (Dt / 4.0) * f.G(M_-1, p_eta_, Q_, beta);
            
            // Step 1.2
            for(int j=0; j<(M_-1); j++)
            {
                coeff = CMt::exp(-(Dt / 8.0) * (p_eta_[j+1] / Q_[j+1]));
                p_eta_[j]*= coeff;
                
                p_eta_[j]+= (Dt / 4.0) * f.G(j, p_eta_, Q_, beta);
                
                coeff = CMt::exp(-(Dt / 8.0) * (p_eta_[j+1] / Q_[j+1]));
                p_eta_[j]*= coeff;
            }
            
            // Step 1.3
            for(int j=0; j<M_; j++)
            {
                eta_[j]+= (Dt / 2.0) * (p_eta_[j] / Q_[j]);
            }
            
            // Step 1.4
            coeff = CMt::exp(-((Dt / 2.0) * (p_eta_[0] / Q_[0])));
            f.scaleMomentum(coeff);
            
            // Step 1.5
            for(int j=(M_-2); j>=0; j--)
            {
                coeff = CMt::exp(-(Dt / 8.0) * (p_eta_[j+1] / Q_[j+1]));
                p_eta_[j]*= coeff;
                
                p_eta_[j]+= (Dt / 4.0) * f.G(j, p_eta_, Q_, beta);
                
                coeff = CMt::exp(-(Dt / 8.0) * (p_eta_[j+1] / Q_[j+1]));
                p_eta_[j]*= coeff;
            }
            
            // Step 1.6
            p_eta_[M_-1]+= (Dt / 4.0) * f.G(M_-1, p_eta_, Q_, beta);
        }
    }
}

void CNHChain::prepareArrays(int N, int dim)
{
    eta_.resize(M_, 0.0);
    p_eta_.resize(M_, 0.0);
    
    Q_.clear();
    Q_.push_back(double(dim*N)*T_*tau_*tau_);
    for(int j=1; j<M_; j++)
        Q_.push_back(T_*tau_*tau_);
}

void CNHChain::setRandNHPos()
{
    for(int j=0; j<M_; j++)
    {
        eta_[j] = double(rand()) / double(RAND_MAX);
    }
}

void CNHChain::setRandNHMom()
{
    double          a = 2.0 / double(RAND_MAX);
    
    for(int j=0; j<M_; j++)
    {
        p_eta_[j] = (a*double(rand()) - 1.0) * (Q_[j] / tau_);
        if(p_eta_[j] == 0.0) p_eta_[j] = (Q_[j] / tau_);
    }
}

END_CUDA_COMPATIBLE()
