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

class CMolTwisterTutorialPool
{
public:
    CMolTwisterTutorialPool() = delete;
    CMolTwisterTutorialPool(std::vector<std::string> possibleTutorialFolderLocations);

public:
    size_t getCount() const;
    std::string getDocumentFilePath(size_t index) const;
    std::string getDocument(size_t index) const;
    std::string getDocumentHeader(size_t index) const;

private:
    std::vector<std::string> documentFileList_;
};