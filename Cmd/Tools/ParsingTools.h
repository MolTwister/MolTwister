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
#include "ToolsBase.h"
#include <string>

class CParsingTools : public CToolsBase
{
public:
    CParsingTools() = delete;
    CParsingTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    bool retrieve1to4BondSepCoeffs(std::vector<std::string> arguments, size_t& arg, double* a1to4BondSepCoeffs) const;
    std::pair<bool, std::string> retrieveDihedralAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, CAtom** atom4Ptr, int* atomIndices) const;
    std::pair<bool, std::string> retrieveAngleAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, int* atomIndices) const;
    std::pair<bool, std::string> retrieveBondAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, int* atomIndices) const;
};
