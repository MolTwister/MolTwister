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
#include "MDFFMatrices.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenCellList
{
public:
    HOSTDEV_CALLABLE CFunctorGenCellList();

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE CMDFFMatrices::CCellListIndex operator()(CMDFFMatrices::CAtom& atom);
    HOST_CALLABLE void assembleCellList(mtdevice_vector<CMDFFMatrices::CCellListIndex>& devAtomCellIndices);
    HOSTDEV_CALLABLE static size_t cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t Nx, size_t Ny);

private:
    int* devCellList_;
    CMDFFMatrices::CListPointer* devCellListEntryPointers_;
    int cellCountX_;
    int cellCountY_;
    int cellCountZ_;
    float pbcWidthX_;
    float pbcWidthY_;
    float pbcWidthZ_;
    float pbcLowX_;
    float pbcLowY_;
    float pbcLowZ_;
};

END_CUDA_COMPATIBLE()
