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
#include <vector>
#include "Utilities/Serializer.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CExpLookup
{
public:
    CExpLookup() { delta_ = 1.0; low_ = 0.0; }
    
public:
    void serialize(CSerializer& io, bool saveToStream);
    void init(double lowArg, double highArg, int N);
    double exp(double arg) const;
    
private:
    double delta_;
    double low_;
    std::vector<double> values_;
};

END_CUDA_COMPATIBLE()
