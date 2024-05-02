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
#include "MolTwisterState.h"
#include "Utilities/ASCIIUtility.h"
#include "Utilities/FileUtility.h"
#include "Utilities/DCDFile.h"
#include "Utilities/3DRect.h"
#include "Utilities/3DBasis.h"
#include "Tools/MolTwisterCmdParser.h"
#include <stdlib.h>
#include <vector>
#include <float.h>

class CCmd
{    
public:
    CCmd() = delete;
    CCmd(CMolTwisterState* state) { commonConstruct(state); parser_ = std::make_shared<CCmdParser>(); }
    virtual ~CCmd() {}
    
public:
    void init() { onRegisterSubCommands(); onAddKeywords(); }
    virtual void execute(std::string commandLine) = 0;
    virtual void onRegisterSubCommands() { }
    bool checkSkipTabRemoveReqFlagAndReset() { bool bRet = skipTabRemoveReqFlag_; skipTabRemoveReqFlag_ = false; return bRet; }
    bool checkCmd(std::string command) const { return (command == getCmd()) ? true : false; }
    virtual std::string getCmd() const = 0;
    virtual std::string getTopLevHelpString() const { return std::string("Help not available for ") + getCmd() + std::string("..."); }
    virtual std::string getHelpString() const { return std::string("\tHelp not available for ") + getCmd() + std::string("!"); }
    virtual std::string getHelpString(std::string subCommand) const { return std::string("\tHelp not available for ") + getCmd() + std::string(" ") + subCommand + std::string("!"); }
    std::shared_ptr<std::vector<std::string>> getListOfSubCommands();
    static int getNumKeywords() { return (int)keywords_.size(); }
    static std::string getKeyword(int index) { return keywords_[index]; }
    void redirectOutput(FILE* stdOut = stdout) { stdOut_ = stdOut; parser_->redirectOutput(stdOut); }
    
protected:
    virtual void onAddKeywords() = 0;
    void addKeyword(std::string keyword);
    void addKeywords(std::vector<std::string> keywords);

private:
    void commonConstruct(CMolTwisterState* state);

protected:
    CMolTwisterState* state_;
    FILE* stdOut_;
    bool skipTabRemoveReqFlag_;
    std::shared_ptr<CCmdParser> parser_;

private:
    static std::vector<std::string> keywords_;
};
