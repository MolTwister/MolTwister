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

class CCmdLoad : public CCmd
{
public:
    CCmdLoad() = delete;
    CCmdLoad(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "load"; }
    std::string getTopLevHelpString() const { return std::string("Load data into MolTwister"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();

private:
    void parseXYZCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseDCDCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseXTCCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parsePDBCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseMTTCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseScriptCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parsePythonCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseMasschargeCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseQeposCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);

    bool readXYZFile(std::string xyzFileName, bool& genBonds);
    bool readDCDFile(std::string dcdFileName, bool& genBonds, int stride);
    bool readXTCFile(std::string xtcFileName, bool& genBonds, int stride);
    bool readPDBFile(std::string pdbFileName, bool& genBonds, const std::pair<bool, std::string>& noQuery);
    bool readMTTFile(std::string mttFileName);
    bool readScriptFile(std::string scriptFileName);
    bool readPythonFile(std::string scriptFileName);
    bool readMassChargeFile(std::string massChargeFileName);
    bool readQEPosFile(std::string qePosFileName, std::string qeInputFileName);
};
