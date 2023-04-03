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
#include <string>

class CFileUtility
{
public:
    CFileUtility() {}
    
public:
    static void removeTabsFromFile(std::string fileName);
    static void addMissingSlash(std::string& folder);
    static std::string getCWD();
    static void removeExtension(std::string& filePath);
    static std::string getExtension(std::string filePath);
    static std::string readLine(FILE* fileHandle, bool& lastLineInFile);
    static std::string readToNextDelimiterIgnoreCppComments(FILE* fileHandle, char delimiter, bool& lastLineInFile);
    static std::string extractLeafDir(std::string absDir);
};
