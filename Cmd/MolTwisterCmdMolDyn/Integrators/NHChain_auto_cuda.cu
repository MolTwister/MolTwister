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
    n_sy = 3;
    n = 4;
    M = 8;
    T = 298.0 / Conv_T;  // [K]
    tau = 20.0 / Conv_t;

    GetSuzukiYoshida(w);
}

void CNHChain::GetSuzukiYoshida(double* w)
{
    const double w0 = 1.0 / (1.0 - pow(2.0, 1.0/3.0));
    
    w[0] = w0;
    w[1] = 1.0 - 2.0*w0;
    w[2] = w0;
}

void CNHChain::Propagator(int N, int dim, double dt, CFct& f)
{
    double  beta = 1.0 / T;
    double  coeff;
    double  pEtaCutoff = 15.0;
    
    
    if((M > 0) && (p_eta[0] > pEtaCutoff))
    {
        printf("Warning! cutoff applied to NH chain (p_eta[0]=%g -> %g)...\r\n", p_eta[0], pEtaCutoff);
        p_eta[0] = 10.0;
    }
    
    for(int mu=0; mu<n_sy; mu++)
    {
        double Dt = w[mu] * dt / double(n);
        for(int i=0; i<n; i++)
        {
            // Step 1.1
            p_eta[M-1]+= (Dt / 4.0) * f.G(M-1, p_eta, Q, beta);
            
            // Step 1.2
            for(int j=0; j<(M-1); j++)
            {
                coeff = CMt::Exp(-(Dt / 8.0) * (p_eta[j+1] / Q[j+1]));
                p_eta[j]*= coeff;
                
                p_eta[j]+= (Dt / 4.0) * f.G(j, p_eta, Q, beta);
                
                coeff = CMt::Exp(-(Dt / 8.0) * (p_eta[j+1] / Q[j+1]));
                p_eta[j]*= coeff;
            }
            
            // Step 1.3
            for(int j=0; j<M; j++)
            {
                eta[j]+= (Dt / 2.0) * (p_eta[j] / Q[j]);
            }
            
            // Step 1.4
            coeff = CMt::Exp(-((Dt / 2.0) * (p_eta[0] / Q[0])));
            f.ScaleMomentum(coeff);
            
            // Step 1.5
            for(int j=(M-2); j>=0; j--)
            {
                coeff = CMt::Exp(-(Dt / 8.0) * (p_eta[j+1] / Q[j+1]));
                p_eta[j]*= coeff;
                
                p_eta[j]+= (Dt / 4.0) * f.G(j, p_eta, Q, beta);
                
                coeff = CMt::Exp(-(Dt / 8.0) * (p_eta[j+1] / Q[j+1]));
                p_eta[j]*= coeff;
            }
            
            // Step 1.6
            p_eta[M-1]+= (Dt / 4.0) * f.G(M-1, p_eta, Q, beta);
        }
    }
}

void CNHChain::PrepareArrays(int N, int dim)
{
    eta.resize(M, 0.0);
    p_eta.resize(M, 0.0);
    
    Q.clear();
    Q.push_back(double(dim*N)*T*tau*tau);
    for(int j=1; j<M; j++)
        Q.push_back(T*tau*tau);
}

void CNHChain::SetRandNHPos()
{
    for(int j=0; j<M; j++)
    {
        eta[j] = double(rand()) / double(RAND_MAX);
    }
}

void CNHChain::SetRandNHMom()
{
    double          a = 2.0 / double(RAND_MAX);
    
    for(int j=0; j<M; j++)
    {
        p_eta[j] = (a*double(rand()) - 1.0) * (Q[j] / tau);
        if(p_eta[j] == 0.0) p_eta[j] = (Q[j] / tau);
    }
}

END_CUDA_COMPATIBLE()
