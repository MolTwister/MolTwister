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
#include "../../Utilities/3DRect.h"
#include "../../Utilities/3DVector.h"
#include "../../Utilities/CUDAGeneralizations.h"
#include <vector>
#include <string>
#include <tuple>

BEGIN_CUDA_COMPATIBLE()

class CMolTwisterStateTools : public CToolsBase
{
public:
    class CCellListEntry
    {
    public:
        CCellListEntry() {}

    public:
        void addAtom(int atomIindex) { atomIndices_.emplace_back(atomIindex); }
        int getNumAtoms() const { return (int)atomIndices_.size(); }
        int getAtom(int index) const { return atomIndices_[index]; }

    private:
        std::vector<int> atomIndices_;
    };

public:
    CMolTwisterStateTools() = delete;
    CMolTwisterStateTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    std::shared_ptr<std::vector<int>> getAtomsWithin(double x, double y, double z, double r, C3DRect pbc, const std::vector<CCellListEntry>& atomCellList, int numPBCBisections) const;
    void generateBonds(double minR, bool verbose=true, bool ignoreHBonds=true, int frame=-1, C3DRect* pbc=nullptr, std::vector<std::string>* bondAtomsToIgnore=nullptr);
    int getNumMolecules() const;
    void getMolecules(std::vector<std::vector<int>>& molecules, const std::vector<std::string>& atomTypesToInclude) const;
    void atomicUnwrap(C3DRect pbc);
    void getAtomsOfMolecule(int molIndex, std::vector<int>& atomIndices) const;
    C3DVector getMoleculeDipoleMoment(const std::vector<int>& atomIndices, int frameIndex, C3DVector& Rc, bool chargeNeutralFormulation) const;
    std::shared_ptr<std::vector<CCellListEntry>> genAtomCellList(C3DRect pbc, int numPBCBisections, int frame) const;
    int xyzToAtomCellIndex(C3DRect pbc, double x, double y, double z, int numPBCBisections) const;
    std::tuple<int, int, int> xyzToAtomCellIndex(C3DRect pbc, double x, double y, double z, int numPBCBisections, bool getCeiling) const;
    void getAllVisibleBondsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2) const;
    void getAllMDBondsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC=false);
    void getAllVisibleAnglesInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3) const;
    void getAllMDAnglesInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC=false);
    void getAllVisibleDihedralsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4) const;
    void getAllMDDihedralsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC=false);
    void reportConsistencyOfMDForceField(FILE* file=stdout);
    C3DVector getCenterOfMass(const std::vector<int>& atomIndices, int frameIndex) const;
    C3DVector getGeometricCenter(const std::vector<int>& atomIndices, int frameIndex) const;
    static C3DVector getGeometricCenter(const std::vector<CAtom*>& atoms, int frameIndex);
    double measureTotCoulombEnergy(const double* a1to4BondSepCoeffs=nullptr, int frame=-1) const;
    double measureCoulombPotential(C3DVector at, const std::vector<std::shared_ptr<CAtom>>* atoms, int frame) const;
    double measureCoulombPotential(C3DVector at, int frame=-1) const;
};

END_CUDA_COMPATIBLE()
