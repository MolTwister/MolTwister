//
// Copyright (C) 2024 Richard Olsen.
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
#include <string>
#include <vector>

class CMolTwisterRestAPI
{
public:
    CMolTwisterRestAPI(const std::vector<std::string>& possibleFolderLocations);

public:
    void run() const;

private:
    static void* threadRun(void* arg);
    static std::string findPythonMainFile(const std::vector<std::string>& possibleFolderLocations);
    static std::string readPythonFile(const std::string& filePath);

private:
    std::string pythonMainFile_;
};
