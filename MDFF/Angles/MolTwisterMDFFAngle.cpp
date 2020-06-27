#include "MolTwisterMDFFAngle.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>
#include <float.h>

BEGIN_CUDA_COMPATIBLE()

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
    double      dTheta_dW = (dDiv == 0.0) ? -1.0 / 1E-10 : -1.0 / dDiv;
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
    double      dTheta_dW = (dDiv == 0.0) ? -1.0 / 1E-10 : -1.0 / dDiv;
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

std::vector<std::pair<float, float>> CMDFFAngle::calc1DForceProfile(float thetaStart, float thetaEnd, int points) const
{
    // We first note that for U=U(theta), the force is F_3=-grad_3 U = -dU/dtheta * (dtheta/dx_3, dtheta/dy_3, dtheta/dz_3),
    // where theta = theta(r_1, r_2, r_3). We which to calculate dU/dtheta. This, is done by calling calcForces()
    // to obtain F_3.x and then calcAngularForceCoeffs13() to obtain dr/dx_3. Thus, (-dU/dtheta) = F_3.x / (dr/dx_3).
    std::vector<std::pair<float, float>> profile;

    if(points <= 0) return profile;

    C3DVector r1(-1.0f, 0.0f, 0.0f);
    C3DVector r2(0.0f, 0.0f, 0.0f);
    C3DVector r3(0.0f, 0.0f, 0.0f);
    C3DVector f1, f2, f3;
    float thetaDelta = (thetaEnd - thetaStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float theta = thetaStart + float(i)*thetaDelta;
        r3 = C3DVector(-cos((double)theta), sin((double)theta), 0.0f);
        calcForces(r1, r2, r3, f1, f2, f3);
        C3DVector c = calcAngularForceCoeffs13(r1, r2, r3);

        profile.emplace_back(std::pair<float, float>(theta, (c.x_ != 0.0f) ? f3.x_ / c.x_ : 0.0f));
    }

    return profile;
}

std::vector<std::pair<float, float>> CMDFFAngle::calc1DPotentialProfile(float thetaStart, float thetaEnd, int points) const
{
    std::vector<std::pair<float, float>> profile;

    if(points <= 0) return profile;

    C3DVector r1(-1.0f, 0.0f, 0.0f);
    C3DVector r2(0.0f, 0.0f, 0.0f);
    C3DVector r3(0.0f, 0.0f, 0.0f);
    float thetaDelta = (thetaEnd - thetaStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float theta = thetaStart + float(i)*thetaDelta;
        r3 = C3DVector(-cos((double)theta), sin((double)theta), 0.0f);
        float E = (float)calcPotential(r1, r2, r3);

        profile.emplace_back(std::pair<float, float>(theta, E));
    }

    return profile;
}

END_CUDA_COMPATIBLE()
