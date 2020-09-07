#pragma once
#include <stdio.h>
#include <vector>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFct
{
public:
    CFct() {}
    virtual ~CFct() {}
    
public:
    virtual double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta) = 0;
    virtual void scaleMomentum(double coeff) = 0;
};

class CNHChain
{
public:
    CNHChain();
    
public:
    void propagator(int N, int dim, double dt, CFct& f);
    void prepareArrays(int N, int dim);
    void setRandNHPos();
    void setRandNHMom();
    
private:
    void getSuzukiYoshida(double* w_);
    
public:
    double T_;                                       // Temperature in reduced units
    int n_;                                          // RESPA steps in Nose-Hoover part
    int M_;                                          // Nose-Hoover chain length
    double tau_;                                     // Thermostat relaxation time in reduced units
    std::vector<double> eta_;                        // Nose-Hoover position coordinates
    std::vector<double> p_eta_;                      // Nose-Hoover momentum
    std::vector<double> Q_;                          // Nose-Hoover Q for each chain entry
    
private:
    int n_sy_;
    double w_[3];
};

END_CUDA_COMPATIBLE()
