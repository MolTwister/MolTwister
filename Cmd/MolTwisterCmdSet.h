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

class CCmdSet : public CCmd
{
public:
    CCmdSet() = delete;
    CCmdSet(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "set"; }
    std::string getTopLevHelpString() const { return std::string("Set various properties of MolTwister and the loaded systems"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();

private:
    void parseProjectionCommand(std::string commandLine, int& arg);
    void parseFullscreenCommand(std::string commandLine, int& arg);
    void parseUserdefpbcCommand(std::string commandLine, int& arg);
    void parseBondacrosspbcCommand(std::string commandLine, int& arg);
    void parseRedrawlimitCommand(std::string commandLine, int& arg);
    void parseFogCommand(std::string commandLine, int& arg);
    void parseUsevdwradiusCommand(std::string commandLine, int& arg);
    void parseVdwscalefactorCommand(std::string commandLine, int& arg);
    void parseLabelsCommand(std::string commandLine, int& arg);
    void parseLabelsfontsizeCommand(std::string commandLine, int& arg);
    void parseBackgroundcolorCommand(std::string commandLine, int& arg);
};
