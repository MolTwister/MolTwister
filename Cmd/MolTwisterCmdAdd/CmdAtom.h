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
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdAtom : public CCmdEntry
{
public:
    CCmdAtom() = delete;
    CCmdAtom(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdAtom() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void addAtomByCubeCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg);
    bool addAtomBySphereCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg);
    C3DVector calcDirVecFromBond(const CAtom* at1, const CAtom* at2, double angleAroundBond, double angleDirBond, double len, int frame) const;
    C3DVector calcDirVecFromAngle(const CAtom* at1, const CAtom* at2, const CAtom* at3,
                                  double angleAroundBond, double angleDirBond, double len, bool& dihedralNotDefined, int frame) const;

private:
    std::string lastError_;
};
