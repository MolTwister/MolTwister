#pragma once
#include <stdio.h>
#include <vector>

using namespace std;

class CFct
{
public:
    CFct() {}
    
public:
    virtual double G(int j, vector<double>& p_eta, vector<double>& Q, double beta) = 0;
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
    vector<double> eta;                             // Nose-Hoover position coordinates
    vector<double> p_eta;                           // Nose-Hoover momentum
    vector<double> Q;                               // Nose-Hoover Q for each chain entry
    
private:
    int n_sy;
    double w[3];
};
