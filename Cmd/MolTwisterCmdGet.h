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

class CCmdGet : public CCmd
{
public:
    CCmdGet() = delete;
    CCmdGet(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "get"; }
    std::string getTopLevHelpString() const { return std::string("Get various properties of MolTwister and the loaded systems"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
    
private:
    void parseAtomtypesCommand(std::string commandLine, int& arg);
    void parseMdinconsistencyCommand(std::string commandLine, int& arg);
    void parseBondinfoCommand(std::string commandLine, int& arg);
    void parseUserdefpbcCommand(std::string commandLine, int& arg);
    void parseDefaultatompropsCommand(std::string commandLine, int& arg);
    void parseGpuinfoCommand(std::string commandLine, int& arg);
};
