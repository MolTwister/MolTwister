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
#include "MolTwisterCmdEntry.h"

class CCmdParser
{
public:
    CCmdParser() = default;

public:
    void registerCmd(std::shared_ptr<CCmdEntry> cmdEntry) { cmdEntryList_.emplace_back(cmdEntry); }
    std::vector<std::string> getCmdLineKeywords();
    std::string executeCmd(std::string commandLine);
    std::string genHelpText(std::string parentCmd);
    std::string genHelpText(std::string parentCmd, std::string subCommand);
    std::shared_ptr<std::vector<std::string>> getListOfCmdEntryCommands();
    void purge() { cmdEntryList_.clear(); }
    void redirectOutput(FILE* stdOut = stdout);

private:
    std::vector<std::shared_ptr<CCmdEntry>> cmdEntryList_;
};
