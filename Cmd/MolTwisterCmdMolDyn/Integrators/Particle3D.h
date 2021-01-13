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

#pragma once
#include "../../../Utilities/3DVector.h"
#include "../../../Utilities/CUDAGeneralizations.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

BEGIN_CUDA_COMPATIBLE()

class CParticle3D
{
public:
    CParticle3D();

public:
    double m_;
    C3DVector r_;
    C3DVector p_;
};

END_CUDA_COMPATIBLE()
