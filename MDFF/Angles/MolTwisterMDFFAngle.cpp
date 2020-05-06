#include "MolTwisterMDFFAngle.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>
#include <float.h>

void CMDFFAngle::parse(std::vector<std::string> arguments)
{    
    size_t nextArg = onParse(arguments);
    comments_ = CASCIIUtility::argsToString(arguments, nextArg);
}

C3DVector CMDFFAngle::calcAngularForceCoeffs13(C3DVector r1, C3DVector r2, C3DVector r3) const
{
    C3DVector   r12 = r1 - r2;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    C3DVector   K = r12*RInv;
    C3DVector   r32 = r3 - r2;
    double      A = K.x_*r32.x_ + K.y_*r32.y_ + K.z_*r32.z_;
    double      B = r32.x_*r32.x_ + r32.y_*r32.y_ + r32.z_*r32.z_;
    double      sqrtB = sqrt(B);
    double      sqrtB3 = (B == 0.0) ? sqrt(1E-30) : sqrtB*sqrtB*sqrtB;
    double      dW_dx = (K.x_*B - r32.x_*A) / sqrtB3;
    double      dW_dy = (K.y_*B - r32.y_*A) / sqrtB3;
    double      dW_dz = (K.z_*B - r32.z_*A) / sqrtB3;
    double      W = (sqrtB == 0.0) ? A / 10E-3 : A / sqrtB;
    double      dDiv = sqrt(1.0 - W*W);
    double      dTheta_dW = (dDiv == 0.0) ? 1.0 / 1E-10 : 1.0 / dDiv;
    C3DVector   dTheta_dr3 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);
    
    return dTheta_dr3;
}

C3DVector CMDFFAngle::calcAngularForceCoeffs2(C3DVector r1, C3DVector r2, C3DVector r3) const
{
    C3DVector   r12 = r1 - r2;
    C3DVector   r32 = r3 - r2;
    double      A = r12.x_*r12.x_ + r12.y_*r12.y_ + r12.z_*r12.z_;
    double      B = r32.x_*r32.x_ + r32.y_*r32.y_ + r32.z_*r32.z_;
    double      C = r12.x_*r32.x_ + r12.y_*r32.y_ + r12.z_*r32.z_;
    double      AB = A*B;
    double      sqrtAB = sqrt(AB);
    double      sqrtAB3 = (AB == 0.0) ? sqrt(1E-30) : sqrtAB*sqrtAB*sqrtAB;
    double      dW_dx = ((2.0*r2.x_ - r1.x_ - r3.x_)*AB + (r32.x_*A + r12.x_*B)*C) / sqrtAB3;
    double      dW_dy = ((2.0*r2.y_ - r1.y_ - r3.y_)*AB + (r32.y_*A + r12.y_*B)*C) / sqrtAB3;
    double      dW_dz = ((2.0*r2.z_ - r1.z_ - r3.z_)*AB + (r32.z_*A + r12.z_*B)*C) / sqrtAB3;
    double      W = (sqrtAB == 0.0) ? A / 10E-3 : C / sqrtAB;
    double      dDiv = sqrt(1.0 - W*W);
    double      dTheta_dW = (dDiv == 0.0) ? 1.0 / 1E-10 : 1.0 / dDiv;
    C3DVector   dTheta_dr2 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);
    
    return dTheta_dr2;
}

bool CMDFFAngle::calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector& f1, C3DVector& f2, C3DVector& f3, double error, int maxIter) const
{
    C3DVector       deltaX, deltaY, deltaZ;
    C3DVector       dU_dr, dU_dr_prev, Diff;
    double          dU, err2 = error*error;
    bool            bFoundReliableResult = true;
    int             i;
    
    
    // Calculate F1
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1+deltaX, r2, r3) - calcPotential(r1, r2, r3);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1+deltaY, r2, r3) - calcPotential(r1, r2, r3);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1+deltaZ, r2, r3) - calcPotential(r1, r2, r3);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f1 = dU_dr;
        Diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(Diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f1*= -1.0;
    if(i >= maxIter) bFoundReliableResult = false;
    
    
    // Calculate F2
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1, r2+deltaX, r3) - calcPotential(r1, r2, r3);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2+deltaY, r3) - calcPotential(r1, r2, r3);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2+deltaZ, r3) - calcPotential(r1, r2, r3);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f2 = dU_dr;
        Diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(Diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f2*= -1.0;
    if(i >= maxIter) bFoundReliableResult = false;

    
    // Calculate F3
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1, r2, r3+deltaX) - calcPotential(r1, r2, r3);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2, r3+deltaY) - calcPotential(r1, r2, r3);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2, r3+deltaZ) - calcPotential(r1, r2, r3);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f3 = dU_dr;
        Diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(Diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f3*= -1.0;
    if(i >= maxIter) bFoundReliableResult = false;
    
    
    return bFoundReliableResult;
}

