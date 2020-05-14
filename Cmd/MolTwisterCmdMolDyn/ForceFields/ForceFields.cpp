#include "ForceFields.h"
#include <math.h>
#include <stdio.h>

CPotLennardJones::CPotLennardJones()
{
    m_bActive = false;
    m_epsilon = 0.0;
    m_sigma6 = 0.0;
    m_sigma12 = 0.0;
}

void CPotLennardJones::Set(double epsilon, double sigma)
{
    double sigma2 = sigma * sigma;
    double sigma3 = sigma2 * sigma;
    m_sigma6 = sigma3 * sigma3;
    m_sigma12 = m_sigma6 * m_sigma6;
    m_epsilon = epsilon;

    m_bActive = true;
}

C3DVector CPotLennardJones::CalcForce(C3DVector r_i, C3DVector r_j)
{
    C3DVector r_vec_ji = r_i - r_j;
    double r_ji = r_vec_ji.norm();
    double r_inv_ji = (r_ji == 0.0) ? 1.0 / 1.0E-3 : 1.0 / r_ji;
    double r_inv_ji2 = r_inv_ji * r_inv_ji;
    double r_inv_ji3 = r_inv_ji2 * r_inv_ji;
    double r_inv_ji6 = r_inv_ji3 * r_inv_ji3;
    double r_inv_ji12 = r_inv_ji6 * r_inv_ji6;
    double c = 4.0*m_epsilon*r_inv_ji2
             * (12.0*m_sigma12*r_inv_ji12 - 6.0*m_sigma6*r_inv_ji6);

    return r_vec_ji * c;
}


void CPotHarmonic::Set(double K, C3DVector r0)
{
    m_K = K;
    m_r0 = r0;

    m_bActive = true;
}

C3DVector CPotHarmonic::CalcForce(C3DVector r_i)
{
    C3DVector r_i0 = r_i - m_r0;
    
    return r_i0 * (-m_K);
}


C3DVector CForceNonBonded::CalcForce(C3DVector r_i, C3DVector r_j)
{
    C3DVector Ftot;
    if(LJ.IsActive()) Ftot+= LJ.CalcForce(r_i, r_j);
    
    return Ftot;
}


C3DVector CForceExternal::CalcForce(C3DVector r_i)
{
    C3DVector Ftot;
    if(Harm.IsActive()) Ftot+= Harm.CalcForce(r_i);
    
    return Ftot;
}


C3DVector CForceHarmBond::CalcForce(C3DVector r_i, C3DVector r_j, double dBoxX, double dBoxY, double dBoxZ)
{
    C3DVector   r_ij = r_i - r_j;
    double      R_ij;
    
    if(fabs(r_ij.x_) > (dBoxX / 2.0))
    {
        if(r_ij.x_ < 0.0) r_ij.x_ = r_ij.x_ + dBoxX;
        else              r_ij.x_ = r_ij.x_ - dBoxX;
    }
    if(fabs(r_ij.y_) > (dBoxY / 2.0))
    {
        if(r_ij.y_ < 0.0) r_ij.y_ = r_ij.y_ + dBoxY;
        else              r_ij.y_ = r_ij.y_ - dBoxY;
    }
    if(fabs(r_ij.z_) > (dBoxZ / 2.0))
    {
        if(r_ij.z_ < 0.0) r_ij.z_ = r_ij.z_ + dBoxZ;
        else              r_ij.z_ = r_ij.z_ - dBoxZ;
    }
    
    R_ij = r_ij.norm();
    
    return r_ij * (-m_K * (R_ij - m_r0) / R_ij);
}

