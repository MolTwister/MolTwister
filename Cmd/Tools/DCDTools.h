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
#include <float.h>
#include "../../Utilities/DCDFile.h"

class CDCDTools : public CToolsBase
{
public:
    class CHBondCriteria
    {
    public:
        CHBondCriteria() { maxLenR_dh_ = DBL_MAX; maxLenR_ha_ = DBL_MAX; maxLenR_da_ = DBL_MAX; minAngle_dha_ = 0.0; }

    public:
        void set(std::string donor, std::string hydrogen, std::string acceptor, double maxLenR_dh, double maxLenR_ha, double maxLenR_da, double minAngle_dha)
            { donor_ = donor; hydrogen_ = hydrogen; acceptor_ = acceptor; maxLenR_dh_ = maxLenR_dh; maxLenR_ha_ = maxLenR_ha; maxLenR_da_ = maxLenR_da; minAngle_dha_ = minAngle_dha; }
        void setDonor(std::string donor) { donor_ = donor; }
        void setHydrogen(std::string hydrogen) { hydrogen_ = hydrogen; }
        void setAcceptor(std::string acceptor) { acceptor_ = acceptor; }
        void setR_dh(double val) { maxLenR_dh_ = val; }
        void setR_ha(double val) { maxLenR_ha_ = val; }
        void setR_da(double val) { maxLenR_da_ = val; }
        void setAlpha_dha(double val) { minAngle_dha_ = val; }
        double R_dh() const { return maxLenR_dh_; }
        double R_ha() const { return maxLenR_ha_; }
        double R_da() const { return maxLenR_da_; }
        double alpha_dha() const { return minAngle_dha_; }
        std::string getDonor() const { return donor_; }
        std::string getHydrogen() const { return hydrogen_; }
        std::string getAcceptor() const { return acceptor_; }
        double calcSizeOfGridToFitHBondPairs() const;

    private:
        // d: donor, h: hydrogen, a: acceptor
        double maxLenR_dh_;
        double maxLenR_ha_;
        double maxLenR_da_;
        double minAngle_dha_;
        std::string donor_;
        std::string hydrogen_;
        std::string acceptor_;
    };

    class CHBond
    {
    public:
        CHBond() { d_ = h_ = a_ = -1; }
        CHBond(int d, int h, int a) { d_ = d; h_ = h; a_ = a; }

    public:
        // d: donor, h: hydrogen, a: acceptor (indices)
        int d_;
        int h_;
        int a_;
    };

    enum EHBondSpan { spanOnlyInterMolecularHBonds=0, spanOnlyIntraMolecularHBonds=1, spanAllHBonds=2 };

public:
    CDCDTools() = delete;
    CDCDTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    void atomicUnwrapCurrDCDFrame(CDCDFile* dcdFile, const std::vector<std::vector<int>>& molList, const C3DRect* pbc) const;
    C3DVector getMoleculeDipoleMoment(const std::vector<int>& atomIndices, const CDCDFile* dcdFile, C3DVector& Rc, bool chargeNeutralFormulation) const;
    C3DVector getCenterOfMass(const std::vector<int>& atomIndices, const CDCDFile* dcdFile) const;
    std::shared_ptr<std::vector<CHBond>> retrieveHBondsCurrDCDFrame(const CDCDFile* dcdFile, const C3DRect* pbc, const std::vector<CHBondCriteria>& HBondCriteria, EHBondSpan HBondSpan, bool noPBC) const;
    static double measureDihedral(const CDCDFile& dcdFile, int index1, int index2, int index3, int index4, const C3DRect* pbc=nullptr);
    static double measureAngle(const CDCDFile& dcdFile, int index1, int index2, int index3, const C3DRect* pbc=nullptr);
    static double measureDist(const CDCDFile& dcdFile, int index1, int index2, const C3DRect* pbc=nullptr);
};
