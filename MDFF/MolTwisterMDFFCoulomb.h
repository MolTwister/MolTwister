#pragma once
#include <stdio.h>
#include "Utilities/3DVector.h"
#include "MolTwisterState.h"

class CMDFFCoulomb
{
public:
    CMDFFCoulomb(CMolTwisterState* molTwisterState) { state_ = molTwisterState; coulombConst_ = 1389.3545752702225; /*[kJAA/mol]*/ }
    
public:
    double calcPotentialBetween(C3DVector r1, C3DVector r2, double q1, double q2) const;
    double calcForceBetween(C3DVector r1, C3DVector r2, double q1, double q2, C3DVector& f1, C3DVector& f2) const;
    bool calcForceBetweenNumerically(C3DVector r1, C3DVector r2, double q1, double q2, C3DVector& f1, C3DVector& f2, double& r12, double error=1E-3, int maxIter=50) const;

private:
    CMolTwisterState* state_;
    
private:
    double coulombConst_;
};
