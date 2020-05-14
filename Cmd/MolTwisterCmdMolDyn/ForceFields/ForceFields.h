#ifndef __ThesisMDTests__Forces__
#define __ThesisMDTests__Forces__

#include "../../../Utilities/3DVector.h"

class CPotLennardJones
{
public:
    CPotLennardJones();
    
public:
    void Set(double epsilon, double sigma);
    bool IsActive() { return m_bActive; }
    C3DVector CalcForce(C3DVector r_i, C3DVector r_j);
    
private:
    bool m_bActive;
    double m_epsilon;
    double m_sigma6;
    double m_sigma12;
};

class CPotHarmonic
{
public:
    CPotHarmonic() { m_bActive = false; m_K = 0.0; }
    
public:
    void Set(double K, C3DVector r0);
    bool IsActive() { return m_bActive; }
    C3DVector CalcForce(C3DVector r_i);
    
public:
    bool m_bActive;
    double m_K;
    C3DVector m_r0;
};

class CForceNonBonded
{
public:
    CForceNonBonded() {}
    
public:
    C3DVector CalcForce(C3DVector r_i, C3DVector r_j);
    
public:
    CPotLennardJones LJ;
};

class CForceExternal
{
public:
    CForceExternal() {}
    
public:
    C3DVector CalcForce(C3DVector r_i);
    
public:
    CPotHarmonic Harm;
};

class CForceHarmBond
{
public:
    CForceHarmBond() { m_K = 0.0; }

public:
    void Set(double K, double r0, int iPart1, int iPart2) { m_K = K; m_i = iPart1; m_j = iPart2; m_r0 = r0; }

public:
    C3DVector CalcForce(C3DVector r_i, C3DVector r_j, double dBoxX, double dBoxY, double dBoxZ);

public:
    double m_K;
    double m_r0;
    int m_i;
    int m_j;
};

#endif
