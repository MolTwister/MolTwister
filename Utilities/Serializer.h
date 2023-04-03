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
#include <vector>
#include <string>
#include <stdio.h>

class CSerializer
{
public:
    CSerializer()
    {
        readPos_ = 0;
    }

public:
    void operator<<(int i)
    {
        archive_.emplace_back(std::to_string(i));
    }

    void operator<<(size_t i)
    {
        archive_.emplace_back(std::to_string(i));
    }

    void operator<<(char c)
    {
        archive_.emplace_back(std::to_string(c));
    }

    void operator<<(double f)
    {
        archive_.emplace_back(std::to_string(f));
    }

    void operator<<(float f)
    {
        archive_.emplace_back(std::to_string(f));
    }

    void operator<<(bool b)
    {
        archive_.emplace_back(std::to_string(b ? 1 : 0));
    }

    void operator<<(std::string s)
    {
        archive_.emplace_back(s);
    }


    void operator>>(int& i)
    {
        i = (int)std::atoi(archive_[readPos_++].data());
    }

    void operator>>(size_t& i)
    {
        i = (size_t)std::atoi(archive_[readPos_++].data());
    }

    void operator>>(char& c)
    {
        c = (char)std::atoi(archive_[readPos_++].data());
    }

    void operator>>(double& f)
    {
        f = (double)std::atof(archive_[readPos_++].data());
    }

    void operator>>(float& f)
    {
        f = (float)std::atof(archive_[readPos_++].data());
    }

    void operator>>(bool& b)
    {
        b = (std::atoi(archive_[readPos_++].data()) == 1) ? true : false;
    }

    void operator>>(std::string& s)
    {
        s = archive_[readPos_++];
    }

    void seekBegin()
    {
        readPos_ = 0;
    }

private:
    std::vector<std::string> archive_;
    int readPos_;
};
