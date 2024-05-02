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

#include "3DBasis.h"

BEGIN_CUDA_COMPATIBLE()

void C3DBasis::generateCartessianBasisAt(C3DVector, C3DVector newBasisZ)
{
    // Calculate Z unit vector of new basis, w, described in old system
    w_ = newBasisZ;
    w_.normalize();

    // Find a vector not parallel to NewBasisZ
    C3DVector l_p, l = newBasisZ, l_pxl;
    l_p = l + C3DVector(1.0, 0.0, 0.0);
    l_pxl = l_p.cross(l);
    if(l_pxl.isZero())
    {
        l_p = l + C3DVector(0.0, 1.0, 0.0);
        l_pxl = l_p.cross(l);
        if(l_pxl.isZero())
        {
            l_p = l + C3DVector(0.0, 0.0, 1.0);
            l_pxl = l_p.cross(l);
        }
    }

    // Calculate X unit vector of new basis, u, described in old system
    u_ = l_pxl;
    u_.normalize();

    // Calculate Y unit vector of new basis, v, described in old system
    v_ = w_.cross(u_);
}

END_CUDA_COMPATIBLE()
