#pragma once
#include <stdio.h>
#include <vector>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFct
{
public:
    CFct() {}
    
public:
    virtual double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta) = 0;
    virtual void ScaleMomentum(double coeff) = 0;
};

class CNHChain
{
public:
    CNHChain();
    
public:
    void Propagator(int N, int dim, double dt, CFct& f);
    void PrepareArrays(int N, int dim);
    void SetRandNHPos();
    void SetRandNHMom();
    
private:
    void GetSuzukiYoshida(double* w);
    
public:
    double T;                                       // Temperature in reduced units
    int n;                                          // RESPA steps in Nose-Hoover part
    int M;                                          // Nose-Hoover chain length
    double tau;                                     // Thermostat relaxation time in reduced units
    std::vector<double> eta;                        // Nose-Hoover position coordinates
    std::vector<double> p_eta;                      // Nose-Hoover momentum
    std::vector<double> Q;                          // Nose-Hoover Q for each chain entry
    
private:
    int n_sy;
    double w[3];
};

END_CUDA_COMPATIBLE()
