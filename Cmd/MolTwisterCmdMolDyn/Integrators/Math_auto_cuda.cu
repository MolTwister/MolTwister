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

#include "Math.h"

#define EXP1 2.71828

BEGIN_CUDA_COMPATIBLE()

double CMt::a_[6] = { 1.0, 1.0/6.0, 1.0/120.0, 1.0/5040.0, 1.0/362880.0, 1.0/39916800.0 };

double CMt::exp(double x)
{
    // Improve numerical stability by linear approximation
    // of exp(x), around x=1, for x>1. Valid for equations
    // where x is small during natural progression.
    if(x > 1.0)
        return EXP1 * x;
    
    return ::exp(x);
}

double CMt::sinhXoverX(double x)
{
    double sum = 0.0;
    double x_n = 1.0;
    
    // Improve numerical stability by rather calculating
    // (1/x)*TaylorExpansion(sinh(x)) to 10th order
    for(int n=0; n<6; n++)
    {
        sum+= a_[n]*x_n;
        x_n*= (x*x);
    }
    
    return sum;
}

END_CUDA_COMPATIBLE()
