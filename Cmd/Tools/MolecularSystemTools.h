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
#include "ToolsBase.h"
#include <vector>
#include <string>
#include <map>
#include "../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMolecularSystemTools : public CToolsBase
{
public:
    CMolecularSystemTools() = delete;
    CMolecularSystemTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    void removeDuplicateBondIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>* auxList=nullptr) const;
    void removeDuplicateAngleIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>* auxList=nullptr) const;
    void removeDuplicateDihIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4, std::vector<int>* auxList=nullptr) const;
    void genMolIndices(const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molIndices) const;
    void getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule) const;

private:
    void getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule, std::map<int, int>& visitationTracker) const;
};

END_CUDA_COMPATIBLE()
