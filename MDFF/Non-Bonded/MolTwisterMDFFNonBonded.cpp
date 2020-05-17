#include "MolTwisterMDFFNonBonded.h"
#include "Utilities/ASCIIUtility.h"
#include <float.h>

bool CMDFFNonBonded::parse(std::vector<std::string> arguments)
{
    std::pair<bool, size_t> nextArg = onParse(arguments);
    if(!nextArg.first) return false;
    comments_ = CASCIIUtility::argsToString(arguments, nextArg.second);

    return true;
}

bool CMDFFNonBonded::calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2, double error, int maxIter) const
{
    C3DVector deltaX, deltaY, deltaZ;
    C3DVector dU_dr, dU_dr_prev, Diff;
    double dU, err2 = error*error;
    bool foundReliableResult = true;
    int i;
    
    
    // Calculate F1
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1+deltaX, r2) - calcPotential(r1, r2);
        dU_dr.x_ = dU / deltaX.x_;

        dU = calcPotential(r1+deltaY, r2) - calcPotential(r1, r2);
        dU_dr.y_ = dU / deltaY.y_;

        dU = calcPotential(r1+deltaZ, r2) - calcPotential(r1, r2);
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
    if(i >= maxIter) foundReliableResult = false;

    
    // Calculate F2
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1, r2+deltaX) - calcPotential(r1, r2);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2+deltaY) - calcPotential(r1, r2);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2+deltaZ) - calcPotential(r1, r2);
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
    if(i >= maxIter) foundReliableResult = false;

    
    return foundReliableResult;
}

C3DVector CMDFFNonBonded::calcNonBondForceCoeffs12(C3DVector r1, C3DVector r2)
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;

    return (r12 * RInv);
}
