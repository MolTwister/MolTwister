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
#include "../Tools/MolTwisterCmdEntry.h"
#include "../../Utilities/FFT1D.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

class CCmdVDOS : public CCmdEntry
{
public:
    CCmdVDOS() = delete;
    CCmdVDOS(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdVDOS() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string fillAtomTrajectoryAvgAtoms(const std::vector<int>& atomIndicesToInclude, int frameFrom, int frameTo, int fftLen, std::vector<std::vector<std::vector<CFFT1D::CCplx>>>& atomTrajectory, CProgressBar& pb, CDCDFile& dcdFile) const;
    std::string fillAtomTrajectoryCOMOfAtoms(const std::vector<int>& atomIndicesToInclude, int frameFrom, int frameTo, int fftLen, std::vector<std::vector<std::vector<CFFT1D::CCplx>>>& atomTrajectory, CProgressBar& pb, CDCDFile& dcdFile) const;

private:
    std::string lastError_;
};
