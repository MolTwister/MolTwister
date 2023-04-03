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
#include <vector>
#include "MolTwisterCmd.h"

class CCmdMeasure : public CCmd
{
public:
    CCmdMeasure() = delete;
    CCmdMeasure(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "measure"; }
    std::string getTopLevHelpString() const { return std::string("Perform a measurement"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;
    
protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
