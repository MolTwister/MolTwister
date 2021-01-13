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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include "ASCIIUtility.h"
#include "FileUtility.h"

void CFileUtility::removeTabsFromFile(std::string fileName)
{
    FILE* file = fopen(fileName.data(), "r");
       
    if(file)
    {
        // Read contents from file
        fseek(file, 0L, SEEK_END);
        long size = ftell(file);
        fseek(file, 0L, SEEK_SET);
        std::string text;
        text.resize(size + 1);
        fread((char*)text.data(), 1, size, file);
        text[size] = '\0';
        if(file) fclose(file);
        
        // Modify contents by removing all '\t' characters
        // from the left margine of the text
        std::string modifiedText;
        modifiedText.resize(size + 1);
        long lIndex = 0;
        long lColumn = 0;
        size = text.size();
        for(long l=0; l<size; l++)
        {
            if((lColumn != 0) || (text[l] != '\t'))
            {
                modifiedText[lIndex++] = text[l];
            }
            
            lColumn++;
            if(text[l] == '\n') lColumn = 0;
        }
        modifiedText[lIndex] = '\0';
        
        // Save contents back to file
        file = fopen(fileName.data(), "w");
        if(file)
        {
            size = modifiedText.size();
            fwrite((char*)modifiedText.data(), 1, size, file);
            fclose(file);
        }
    }
}

void CFileUtility::addMissingSlash(std::string& folder)
{
    int strLen = (int)folder.size();
    int i = strLen - 1;
    char character;

    do
    {
        character = folder[i--];

    } while((character == ' ') || (character == '\0'));

    if(character == '/') return;

    folder += "/";
}

std::string CFileUtility::getCWD()
{
    std::string cwd;

    char* str = new char[4096];
    getcwd(str, 4095);
    cwd = str;
    if(str) delete [] str;

    addMissingSlash(cwd);

    return cwd;
}

void CFileUtility::removeExtension(std::string& filePath)
{
    if(filePath.size() == 0) return;
    size_t dotIndex = filePath.rfind('.');
    if(dotIndex == std::string::npos) return;
    filePath = filePath.substr(0, dotIndex);
}

std::string CFileUtility::getExtension(std::string filePath)
{
    if(filePath.size() == 0) return "";
    size_t dotIndex = filePath.rfind('.');
    if(dotIndex == std::string::npos) return "";
    size_t afterDotIndex = dotIndex + 1;
    if(afterDotIndex >= filePath.size()) return "";
    return filePath.substr(afterDotIndex, std::string::npos);
}

std::string CFileUtility::readLine(FILE* fileHandle, bool& lastLineInFile)
{
    char character;
    size_t readCnt;

    if(!fileHandle) return "";

    std::string line;
    do
    {
        readCnt = fread(&character, 1, 1, fileHandle);
        if((character != '\r') && (character != '\n'))
            line+= character;

    } while((readCnt > 0) && (character != '\n'));

    if(readCnt < 1) lastLineInFile = true;
    else            lastLineInFile = false;

    // It could be that this was still the last complete line, but
    // that it ended with <line feed>. So, check next character to see if
    // it is the last in the file
    if(!lastLineInFile)
    {
        fpos_t lastPos;
        fgetpos(fileHandle, &lastPos);
        readCnt = fread(&character, 1, 1, fileHandle);
        if(readCnt < 1) lastLineInFile = true;
        fsetpos(fileHandle, &lastPos);
    }

    return line;
}

std::string CFileUtility::readToNextDelimiterIgnoreCppComments(FILE* fileHandle, char delimiter, bool& lastLineInFile)
{
    char character;
    long ptToFilePos;
    size_t readCnt;

    if(!fileHandle) return "";
    std::string line;

    do
    {
        readCnt = fread(&character, 1, 1, fileHandle);
        if(character == '/')
        {
            ptToFilePos = ftell(fileHandle);
            readCnt = fread(&character, 1, 1, fileHandle);
            if((readCnt > 0) && (character == '/'))
            {
                do
                {
                    readCnt = fread(&character, 1, 1, fileHandle);

                } while((readCnt > 0) && (character != '\n'));

                readCnt = fread(&character, 1, 1, fileHandle);
            }
            else if((readCnt > 0) && (character == '*'))
            {
                char prevChar = ' ';
                bool endComment = false;
                do
                {
                    readCnt = fread(&character, 1, 1, fileHandle);
                    if((character == '/') && (prevChar == '*')) endComment = true;
                    prevChar = character;

                } while((readCnt > 0) && !endComment);

                readCnt = fread(&character, 1, 1, fileHandle);
            }
            else
            {
                fseek(fileHandle, ptToFilePos, SEEK_SET);
            }
        }
        if((character != delimiter) && (character != '\r') && (character != '\n'))
            line+= character;

    } while((readCnt > 0) && (character != delimiter));

    if(readCnt < 1) lastLineInFile = true;
    else            lastLineInFile = false;

    return line;
}

std::string CFileUtility::extractLeafDir(std::string absDir)
{
    std::string szAbsDirCopy = absDir;
    int iStrLen, iSlashIndex = 0;

    CASCIIUtility::removeWhiteSpace(szAbsDirCopy);

    iStrLen = (int)szAbsDirCopy.length();
    if(szAbsDirCopy[iStrLen-1] == '/') szAbsDirCopy.erase(iStrLen-1, 1);

    iStrLen = (int)szAbsDirCopy.length();
    for(int i=iStrLen-1; i>=0; i--)
    {
        if(szAbsDirCopy[i] == '/') break;
        iSlashIndex = i;
    }

    return szAbsDirCopy.substr(iSlashIndex, std::string::npos);
}
