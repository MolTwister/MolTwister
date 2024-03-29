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

class CCmdPrint : public CCmd
{
public:
    CCmdPrint() = delete;
    CCmdPrint(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "print"; }
    std::string getTopLevHelpString() const { return std::string("Print information (e.g., atom pos. in different file formats)"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseXYZCommand(std::string commandLine, int& arg);
    void parsePDBCommand(std::string commandLine, int& arg);
    void parseMTTCommand(std::string commandLine, int& arg);
    void parseVersionCommand(std::string commandLine, int& arg);
    void parseMixffCommand(std::string commandLine, int& arg);
    bool getMixFFConstituents(std::string commandLine, int& arg, std::vector<std::string>& constituents) const;
};
