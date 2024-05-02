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

#include <iostream>
#include <math.h>
#include "ExpLookup.h"

BEGIN_CUDA_COMPATIBLE()

void CExpLookup::serialize(CSerializer& io, bool saveToStream)
{
    if(saveToStream)
    {
        io << values_.size();
        for(double val : values_)
        {
            io << val;
        }

        io << delta_;
        io << low_;
    }
    else
    {
        size_t size;
        io >> size;
        values_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            double val;
            io >> val;
            values_[i] = val;
        }

        io >> delta_;
        io >> low_;
    }
}

void CExpLookup::init(double lowArg, double highArg, int N)
{
    if(N == 0) return;
    
    values_.reserve(N);
    
    delta_ = (highArg - lowArg) / double(N);
    low_ = lowArg;
    
    for(int i=0; i<N; i++)
    {
        values_.emplace_back(::exp(lowArg + delta_ * double(i)));
    }
}

double CExpLookup::exp(double arg) const
{
    int n = int((arg - low_) / delta_);

    if((n < 0) || (n >= (int)values_.size())) return ::exp(arg);
    
    return values_[n];
}

END_CUDA_COMPATIBLE()
