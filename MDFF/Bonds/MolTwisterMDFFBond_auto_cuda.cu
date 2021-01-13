//
// Copyright (C) 2021 Richard Olsen.
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

#include "MolTwisterMDFFBond.h"
#include "Utilities/ASCIIUtility.h"
#include <float.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFBond::serialize(CSerializer& io, bool saveToStream)
{
    int bondDetectionCriteria;

    if(saveToStream)
    {
        bondDetectionCriteria = (int)bondDetectionCriteria_;

        io << atomNamesToBond_[0];
        io << atomNamesToBond_[1];
        io << comments_;
        io << bondDetectionCriteria;
        io << detCritR0_;
    }
    else
    {
        io >> atomNamesToBond_[0];
        io >> atomNamesToBond_[1];
        io >> comments_;
        io >> bondDetectionCriteria;
        io >> detCritR0_;

        bondDetectionCriteria_ = (EDetCrit)bondDetectionCriteria;
    }
}

void CMDFFBond::parse(std::vector<std::string> arguments)
{
    size_t nextArg = onParse(arguments);
    comments_ = CASCIIUtility::argsToString(arguments, nextArg);
}

bool CMDFFBond::calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2, double error, int maxIter) const
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
    for(i=0; i<maxIter; i++)
    {
        dU = calcPotential(r1+deltaX, r2) - calcPotential(r1, r2);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1+deltaY, r2) - calcPotential(r1, r2);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1+deltaZ, r2) - calcPotential(r1, r2);
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
        dU = calcPotential(r1, r2+deltaX) - calcPotential(r1, r2);
        dU_dr.x_ = dU / deltaX.x_;
        
        dU = calcPotential(r1, r2+deltaY) - calcPotential(r1, r2);
        dU_dr.y_ = dU / deltaY.y_;
        
        dU = calcPotential(r1, r2+deltaZ) - calcPotential(r1, r2);
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

std::vector<std::pair<float, float>> CMDFFBond::calc1DForceProfile(float rStart, float rEnd, int points) const
{
    // We first note that for U=U(r), the force is F_2=-grad_2 U = -dU/dr * (dr/dx_2, dr/dy_2, dr/dz_2),
    // where r = sqrt((x_2 - x_1)^2 + (y_2 - y_1)^2 + (z_2 - z_1)^2). We which to calculate (-dU/dr). This,
    // is done by calling calcForces() to obtain F_2.x and then calcBondForceCoeffs12() to obtain dr/dx_2.
    // Thus, (-dU/dr) = F_2.x / (dr/dx_2).
    std::vector<std::pair<float, float>> profile;

    if(points <= 0) return profile;

    C3DVector r1(0.0, 0.0, 0.0);
    C3DVector r2(0.0, 0.0, 0.0);
    C3DVector f1, f2;
    float rDelta = (rEnd - rStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float r = rStart + float(i)*rDelta;
        r2.x_ = double(r);
        calcForces(r1, r2, f1, f2);
        C3DVector c = calcBondForceCoeffs12(r1, r2);

        profile.emplace_back(std::pair<float, float>(r, (c.x_ != 0.0) ? f2.x_ / c.x_ : 0.0));
    }

    return profile;
}

std::vector<std::pair<float, float>> CMDFFBond::calc1DPotentialProfile(float rStart, float rEnd, int points) const
{
    std::vector<std::pair<float, float>> profile;

    if(points <= 0) return profile;

    C3DVector r1(0.0, 0.0, 0.0);
    C3DVector r2(0.0, 0.0, 0.0);
    float rDelta = (rEnd - rStart) / float(points-1);
    for(int i=0; i<points; i++)
    {
        float r = rStart + float(i)*rDelta;
        r2.x_ = double(r);
        float E = (float)calcPotential(r1, r2);

        profile.emplace_back(std::pair<float, float>(r, E));
    }

    return profile;
}

C3DVector CMDFFBond::calcBondForceCoeffs12(C3DVector r1, C3DVector r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;

    return (r12 * RInv);
}

END_CUDA_COMPATIBLE()
