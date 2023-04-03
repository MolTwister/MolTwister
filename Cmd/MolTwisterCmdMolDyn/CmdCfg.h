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
#include "Config/MolDynConfig.h"

class CCmdCfg : public CCmdEntry
{
public:
    CCmdCfg() = delete;
    CCmdCfg(CMolTwisterState* state, FILE* stdOut, CMolDynConfig* molDynConfig) : CCmdEntry(state, stdOut) { molDynConfig_ = molDynConfig; }
    virtual ~CCmdCfg() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void parseGetCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseSetCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseResettodefaultsCommand(const std::vector<std::string>& arguments, size_t& arg);

private:
    std::string lastError_;
    CMolDynConfig* molDynConfig_;
};
