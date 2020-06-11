#include "MolTwisterMDFFDih.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>
#include <float.h>

void CMDFFDih::parse(std::vector<std::string> arguments)
{
    size_t nextArg = onParse(arguments);
    comments_ = CASCIIUtility::argsToString(arguments, nextArg);
}

C3DVector CMDFFDih::calcDihedralForceCoeffs14(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    C3DVector       n1 = r23.cross(r21);
    double          n1Norm = n1.norm();
    double          n1NormInv = (n1Norm == 0.0) ? 1.0 / 1E-10 : 1.0 / n1Norm;
    C3DVector       K = n1 * n1NormInv;
    double          a = r32.y_*r34.z_ - r32.z_*r34.y_;
    double          b = r32.z_*r34.x_ - r32.x_*r34.z_;
    double          c = r32.x_*r34.y_ - r32.y_*r34.x_;
    double          A = K.x_*a + K.y_*b + K.z_*c;
    double          B = a*a + b*b + c*c;
    double          sqrtBInv = (B == 0.0) ? 1.0 / 1E-20 : 1.0 / sqrt(B);
    double          sqrtBInv3 = sqrtBInv*sqrtBInv*sqrtBInv;
    double          W = A*sqrtBInv;
    double          dW_dx = ((K.y_*r32.z_ - K.z_*r32.y_)*B - (r32.z_*(r32.z_*r34.x_ - r32.x_*r34.z_) - r32.y_*(r32.x_*r34.y_ - r32.y_*r34.x_))*A) * sqrtBInv3;
    double          dW_dy = ((K.z_*r32.x_ - K.x_*r32.z_)*B - (r32.x_*(r32.x_*r34.y_ - r32.y_*r34.x_) - r32.z_*(r32.y_*r34.z_ - r32.z_*r34.y_))*A) * sqrtBInv3;
    double          dW_dz = ((K.x_*r32.y_ - K.y_*r32.x_)*B - (r32.y_*(r32.y_*r34.z_ - r32.z_*r34.y_) - r32.x_*(r32.z_*r34.x_ - r32.x_*r34.z_))*A) * sqrtBInv3;
    double          div = sqrt(1.0 - W*W);
    double          dTheta_dW = (div == 0.0) ? 1.0 / 1E-10 : 1.0 / div;
    C3DVector       dTheta_dr4 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);
    
    return dTheta_dr4;
}

C3DVector CMDFFDih::calcDihedralForceCoeffs23(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    double          K1 = r32.z_*r34.x_ - r32.x_*r34.z_;
    double          K2 = r34.z_ - r32.z_;
    double          K3 = r23.z_*r21.x_ - r23.x_*r21.z_;
    double          K4 = r32.x_*r34.y_ - r32.y_*r34.x_;
    double          K5 = r32.y_ - r34.y_;
    double          K6 = r23.x_*r21.y_ - r23.y_*r21.x_;
    double          L1 = r32.y_*r34.z_ - r32.z_*r34.y_;
    double          L2 = -K2;
    double          L3 = r23.y_*r21.z_ - r23.z_*r21.y_;
    double          L4 = K4;
    double          L5 = r34.x_ - r32.x_;
    double          L6 = K6;
    double          M1 = L1;
    double          M2 = -K5;
    double          M3 = L3;
    double          M4 = K1;
    double          M5 = -L5;
    double          M6 = K3;
    double          A = L3*L3 + K3*K3 + K6*K6;
    double          B = L1*L1 + K1*K1 + K4*K4;
    double          C = L3*L1 + K3*K1 + K6*K4;
    double          AB = A*B;
    double          sqrtABInv = (AB == 0.0) ? 1.0 / 1E-20 : 1.0 / sqrt(AB);
    double          sqrtABInv3 = sqrtABInv*sqrtABInv*sqrtABInv;
    double          W = C * sqrtABInv;
    double          dC_dx = K2*K3 + r21.y_*K4 + K5*K6 - r21.z_*K1;
    double          dC_dy = L2*L3 - r21.x_*L4 + L5*L6 + r21.z_*L1;
    double          dC_dz = M2*M3 + r21.x_*M4 + M5*M6 - r21.y_*M1;
    double          dW_dx = (dC_dx*AB - C*((r21.y_*K6 - r21.z_*K3)*B + (K1*K2 + K4*K5)*A)) * sqrtABInv3;
    double          dW_dy = (dC_dy*AB - C*((r21.z_*L3 - r21.x_*L6)*B + (L1*L2 + L4*L5)*A)) * sqrtABInv3;
    double          dW_dz = (dC_dz*AB - C*((r21.x_*M6 - r21.y_*M3)*B + (M1*M2 + M4*M5)*A)) * sqrtABInv3;
    double          div = sqrt(1.0 - W*W);
    double          dTheta_dW = (div == 0.0) ? 1.0 / 1E-10 : 1.0 / div;
    C3DVector       dTheta_dr3 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);
    
    return dTheta_dr3;
}

bool CMDFFDih::calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4, C3DVector& f1, C3DVector& f2, C3DVector& f3, C3DVector& f4, double error, int iMaxIter) const
{
    C3DVector deltaX, deltaY, deltaZ;
    C3DVector dU_dr, dU_dr_prev, diff;
    double dU, err2 = error*error;
    bool foundReliableResult = true;
    int i;
    

    // Calculate F1
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<iMaxIter; i++)
    {
        dU = calcPotential(r1+deltaX, r2, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1+deltaY, r2, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1+deltaZ, r2, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f1 = dU_dr;
        diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f1*= -1.0;
    if(i >= iMaxIter) foundReliableResult = false;
    
    
    // Calculate F2
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<iMaxIter; i++)
    {
        dU = calcPotential(r1, r2+deltaX, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2+deltaY, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2+deltaZ, r3, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f2 = dU_dr;
        diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f2*= -1.0;
    if(i >= iMaxIter) foundReliableResult = false;
    
    
    // Calculate F3
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<iMaxIter; i++)
    {
        dU = calcPotential(r1, r2, r3+deltaX, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2, r3+deltaY, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2, r3+deltaZ, r4) - calcPotential(r1, r2, r3, r4);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f3 = dU_dr;
        diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f3*= -1.0;
    if(i >= iMaxIter) foundReliableResult = false;

    
    // Calculate F4
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<iMaxIter; i++)
    {
        dU = calcPotential(r1, r2, r3, r4+deltaX) - calcPotential(r1, r2, r3, r4);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2, r3, r4+deltaY) - calcPotential(r1, r2, r3, r4);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2, r3, r4+deltaZ) - calcPotential(r1, r2, r3, r4);
        dU_dr.z_ = dU / deltaZ.z_;
        
        f4 = dU_dr;
        diff = dU_dr - dU_dr_prev;
        dU_dr_prev = dU_dr;
        if(diff.norm2() < err2)
            break;
        
        deltaX*= 0.5;
        deltaY*= 0.5;
        deltaZ*= 0.5;
    }
    
    f4*= -1.0;
    if(i >= iMaxIter) foundReliableResult = false;
    
    
    return foundReliableResult;
}

std::vector<std::pair<float, float>> CMDFFDih::calc1DForceProfile(float phiStart, float phiEnd, int points) const
{
    // We first note that for U=U(phi), the force is F_4=-grad_4 U = -dU/dphi * (dphi/dx_4, dphi/dy_4, dphi/dz_4),
    // where phi = phi(r_1, r_2, r_3, r4). We which to calculate dU/dphi. This, is done by calling calcForces()
    // to obtain F_4.x and then calcDihedralForceCoeffs14() to obtain dr/dx_4. Thus, (-dU/dphi) = F_4.x / (dr/dx_4).
    std::vector<std::pair<float, float>> profile;

    if(points <= 0) return profile;

    C3DVector r1(0.0f, 1.0f, 0.0f);
    C3DVector r2(0.0f, 0.0f, 0.0f);
    C3DVector r3(0.0f, 0.0f, 1.0f);
    C3DVector r4(0.0f, 1.0f, 1.0f);
    C3DVector f1, f2, f3, f4;
    float phiDelta = (phiEnd - phiStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float phi = phiStart + float(i)*phiDelta;
        r4 = C3DVector(sin((double)phi), cos((double)phi), 1.0f);
        calcForces(r1, r2, r3, r4, f1, f2, f3, f4);
        C3DVector c = calcDihedralForceCoeffs14(r1, r2, r3, r4);

        profile.emplace_back(std::pair<float, float>(phi, (c.x_ != 0.0f) ? f4.x_ / c.x_ : 0.0f));
    }

    return profile;
}

std::vector<std::pair<float, float>> CMDFFDih::calc1DPotentialProfile(float phiStart, float phiEnd, int points) const
{
    std::vector<std::pair<float, float>> profile;

    C3DVector r1(0.0f, 1.0f, 0.0f);
    C3DVector r2(0.0f, 0.0f, 0.0f);
    C3DVector r3(0.0f, 0.0f, 1.0f);
    C3DVector r4(0.0f, 1.0f, 1.0f);
    C3DVector f1, f2, f3, f4;
    float phiDelta = (phiEnd - phiStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float phi = phiStart + float(i)*phiDelta;
        r4 = C3DVector(sin((double)phi), cos((double)phi), 1.0f);
        float E = (float)calcPotential(r1, r2, r3, r4);

        profile.emplace_back(std::pair<float, float>(phi, E));
    }

    return profile;
}
