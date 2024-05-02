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
#include "3DVector.h"
#include "CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class C3DBasis
{
public:
    C3DBasis() {}
    C3DBasis(C3DVector u, C3DVector v, C3DVector w) { u_ = u; v_ = v; w_ = w; }

public:
    void generateCartessianBasisAt(C3DVector pos, C3DVector newBasisZ);

public:
    C3DVector u_;
    C3DVector v_;
    C3DVector w_;
};

END_CUDA_COMPATIBLE()
