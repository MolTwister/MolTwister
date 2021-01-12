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
#include "../../../Utilities/CUDAGeneralizations.h"
#include "MDFFMatrices.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenNeighList
{
public:
    HOSTDEV_CALLABLE CFunctorGenNeighList();

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE int operator()(CMDFFMatrices::CAtom& atom);
    HOSTDEV_CALLABLE static size_t neighIndexToFlatIndex(size_t atomIndex, size_t neighIndex, int maxNeighbors);
    int checkRequiredMaxNumNeighbours(mtdevice_vector<int>& devNeighListCount);

private:
    int* devCellList_;
    CMDFFMatrices::CListPointer* devCellListEntryPointers_;
    CMDFFMatrices::CCellListIndex* devAtomCellIndicesRaw_;
    int* devNeighList_;
    int* devNeighListCount_;
    CMDFFMatrices::CAtom* devAtomList_;
    int numAtoms_;
    int cellCountX_;
    int cellCountY_;
    int cellCountZ_;
    int maxNeighbors_;
    float rCutoff2_;
    C3DRect pbc_;
};

END_CUDA_COMPATIBLE()
