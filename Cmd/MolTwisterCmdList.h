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
#include "MolTwisterCmd.h"

class CCmdList : public CCmd
{
public:
    CCmdList() = delete;
    CCmdList(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "list"; }
    std::string getTopLevHelpString() const { return std::string("List atomic properties"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();    

private:
    void parseAllCommand(std::string commandLine, int& arg);
    void parseMolCommand(std::string commandLine, int& arg);
    void parseFFCommand(std::string commandLine, int& arg);
    void parseLatexCommand(std::string commandLine, int& arg);

    void printListAtomicItem(const CAtom* atom) const;
    void printListNonBondedItem(const CMDFFNonBonded* nonBonded, int index) const;
    void printListBondItem(const CMDFFBond* bond, int index) const;
    void printListAngleItem(const CMDFFAngle* angle, int index) const;
    void printListDihItem(const CMDFFDih* dih, int index) const;
    void printListHeaderAtoms() const;
    void printListHeaderMDNonBond() const;
    void printListHeaderMDBonds() const;
    void printListHeaderMDAngles() const;
    void printListHeaderMDDihedrals() const;

    void printLaTeXListCharges(bool useLongTable) const;
    
    void printLaTeXListNonBondedItem(const CMDFFNonBonded* nonBonded, int index, bool useLongTable) const;
    void printLaTeXListBondItem(const CMDFFBond* bond, int index, bool useLongTable) const;
    void printLaTeXListAngleItem(const CMDFFAngle* angle, int index, bool useLongTable) const;
    void printLaTeXListDihItem(const CMDFFDih* dih, int index, bool useLongTable) const;
    void printLaTeXListHeaderMDNonBond(bool useLongTable) const;
    void printLaTeXListFooterMDNonBond(bool useLongTable) const;
    void printLaTeXListHeaderMDBonds(bool useLongTable) const;
    void printLaTeXListFooterMDBonds(bool useLongTable) const;
    void printLaTeXListHeaderMDAngles(bool useLongTable) const;
    void printLaTeXListFooterMDAngles(bool useLongTable) const;
    void printLaTeXListHeaderMDDihedrals(bool useLongTable) const;
    void printLaTeXListFooterMDDihedrals(bool useLongTable) const;
};
