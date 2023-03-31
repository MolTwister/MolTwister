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

#include <iostream>
#include <filesystem>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "MolTwisterTutorialPool.h"

CMolTwisterTutorialPool::CMolTwisterTutorialPool(std::vector<std::string> possibleTutorialFolderLocations)
{
    for(const std::string& tutorialFolder : possibleTutorialFolderLocations)
    {
        struct stat info;
        if((stat(tutorialFolder.data(), &info) == 0) && (info.st_mode & S_IFDIR))
        {
            fflush(stdout);
            for(const auto& entry : std::filesystem::directory_iterator(tutorialFolder))
            {
                documentFileList_.emplace_back(std::filesystem::absolute(entry.path()));
            }
            break;
        }
    }
}

size_t CMolTwisterTutorialPool::getCount() const
{
    return documentFileList_.size();
}

std::string CMolTwisterTutorialPool::getDocumentFilePath(size_t index) const
{
    return documentFileList_[index];
}

std::string CMolTwisterTutorialPool::getDocument(size_t index) const
{
    std::string document;
    FILE* file = fopen(documentFileList_[index].c_str(), "r");
    if(file)
    {
        fseek(file , 0 , SEEK_END);
        long size = ftell(file);
        rewind(file);

        char* buffer = new char[size];
        fread(buffer, 1, size, file);
        document = buffer;
        if(buffer) delete [] buffer;

        fclose(file);
    }

    return document;
}

std::string CMolTwisterTutorialPool::getDocumentHeader(size_t index) const
{
    std::string document = getDocument(index);
    std::string documentFilePath = getDocumentFilePath(index);
    const std::string emptyMsg = std::string("Tutorial document ") + documentFilePath + std::string(" is empty!");

    if(document.size() == 0) return emptyMsg;
    int startIndexOfMainHeader = document.find("# ");
    if(startIndexOfMainHeader == std::string::npos) return "N/A";
    startIndexOfMainHeader+= 2;

    int endIndexOfMainHeader = document.find("\n", startIndexOfMainHeader + 1);
    if(endIndexOfMainHeader == std::string::npos) endIndexOfMainHeader = document.size() - 1;
    endIndexOfMainHeader--;

    int headerLength = int(endIndexOfMainHeader) - int(startIndexOfMainHeader) + 1;
    if(headerLength < 0) return emptyMsg;
    return document.substr(startIndexOfMainHeader, headerLength);
}
