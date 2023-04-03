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
#include "../../MolTwisterAtom.h"
#include "../../Utilities/3DRect.h"
#include "../../Utilities/3DBasis.h"
#include "../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMolecularTools : public CToolsBase
{
public:
    CMolecularTools() = delete;
    CMolecularTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    static double measureDihedral(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, const CAtom* atom4, int frame, const C3DRect* pbc=nullptr);
    static double measureAngle(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, int frame, const C3DRect* pbc=nullptr);
    static double measureDist(const CAtom* atom1, const CAtom* atom2, int frame, const C3DRect* pbc=nullptr);
    double calcDistance2(const CAtom* fromAtom, const CAtom* toAtom, int frame, const C3DRect& pbc, bool distAcrossPBC) const;
    static C3DVector rotatePosAroundBasisW(C3DVector basisPos, C3DBasis basis, C3DVector posToRotate, double relAngleToRotate);
    static void modBondLengthTo(const CAtom* atom1, CAtom* atom2, double dist, int frame);
    static void modBondTypeTo(CAtom *atom1, CAtom* atom2, const std::string& type, int frame);
    static void modAngleTo(const CAtom* atom1, const CAtom* atom2, CAtom* atom3, double angle, int frame);
    static void modDihedralTo(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, CAtom* atom4, double angle, int frame);
};

END_CUDA_COMPATIBLE()
