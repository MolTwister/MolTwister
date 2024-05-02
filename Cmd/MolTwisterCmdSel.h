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
#include "MolTwisterCmd.h"

class CCmdSel : public CCmd
{
public:
    CCmdSel() = delete;
    CCmdSel(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "sel"; }
    std::string getTopLevHelpString() const { return std::string("Create a selection of atoms"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseAtomCommand(std::string commandLine, int& arg);
    void parseAtomnameCommand(std::string commandLine, int& arg);
    void parseAllCommand(std::string commandLine, int& arg);
    void parseNoneCommand(std::string commandLine, int& arg);
};
