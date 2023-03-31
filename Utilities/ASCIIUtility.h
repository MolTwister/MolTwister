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
#include <string>
#include <vector>
#include <float.h>
#include <climits>

class CASCIIUtility
{
public:
    static std::string extractString(int startIndex, int endIndex, std::string stringIn);
    static bool isWhiteSpace(char character, const char* whiteSpaceChars=nullptr);
    static void removeWhiteSpace(std::string& string, const char* whiteSpaceChars=nullptr);
    static std::string trimFromLeft(const std::string& str, const char* whiteSpaceChars=nullptr, int maxTrim=INT_MAX);
    static std::string trimFromRight(const std::string& str, const char* whiteSpaceChars=nullptr, int maxTrim=INT_MAX);
    static std::string getWord(std::string line, int wordIndex, const char* whiteSpaceChars=nullptr);
    static std::vector<std::string> getWords(std::string line, const char* whiteSpaceChars=nullptr);
    static std::vector<std::string> getLines(std::string text);
    static std::string getDelimitedWord(std::string line, int wordIndex, char startDelimiter, char endDelimiter);
    static std::string removeCRLF(std::string text);
    static int findString(std::string string, std::string line);
    static void parseCLikeFuncCall(std::string inputString, std::string& funcName, std::string& arguments);
    static bool isLineEmpty(std::string line);
    static std::string getArg(const std::vector<std::string>& arguments, size_t argIndex);
    static std::string argsToString(const std::vector<std::string>& arguments);
    static std::string argsToString(const std::vector<std::string>& arguments, size_t firstArgToInclude);
    static std::string createMarkDownCodeBlock(std::string str, int numSpaces, bool removeFirstTab=false);
    static std::string addTabsToDocument(const std::string& document);
};
