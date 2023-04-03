//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include "MolTwisterMDFFCoulomb.h"
#include <float.h>

double CMDFFCoulomb::calcPotentialBetween(C3DVector r1, C3DVector r2, double q1, double q2) const
{
    C3DVector   r12 = r2 - r1;
    double      R12 = r12.norm();
    double      R12Inv = (R12 == 0.0) ? 1.0 / 1E-10 : 1.0 / R12;
    
    return coulombConst_ * q1*q2 * R12Inv;
}

double CMDFFCoulomb::calcForceBetween(C3DVector r1, C3DVector r2, double q1, double q2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv_dUdr_neg = coulombConst_ * q1*q2 * RInv*RInv*RInv;
    C3DVector   F = r12*RInv_dUdr_neg;
    
    f2 = F;
    f1 = F*(-1.0);
    
    return R;
}

bool CMDFFCoulomb::calcForceBetweenNumerically(C3DVector r1, C3DVector r2, double q1, double q2, C3DVector& f1, C3DVector& f2, double&, double error, int maxIter) const
{
    C3DVector       deltaX, deltaY, deltaZ;
    C3DVector       dU_dr, dU_dr_prev, diff;
    double          dU, err2 = error*error;
    bool            foundReliableResult = true;
    int             i;    
    

    // Calculate F1
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotentialBetween(r1+deltaX, r2, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotentialBetween(r1+deltaY, r2, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotentialBetween(r1+deltaZ, r2, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
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
    if(i >= maxIter) foundReliableResult = false;
    
    
    // Calculate F2
    dU_dr_prev = C3DVector(DBL_MAX, DBL_MAX, DBL_MAX);
    deltaX = C3DVector(0.02, 0.0, 0.0);
    deltaY = C3DVector(0.0, 0.02, 0.0);
    deltaZ = C3DVector(0.0, 0.0, 0.02);
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotentialBetween(r1, r2+deltaX, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotentialBetween(r1, r2+deltaY, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotentialBetween(r1, r2+deltaZ, q1, q2) - calcPotentialBetween(r1, r2, q1, q2);
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
    if(i >= maxIter) foundReliableResult = false;


    return foundReliableResult;
}
