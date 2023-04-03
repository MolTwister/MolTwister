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
