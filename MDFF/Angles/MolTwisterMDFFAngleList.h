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
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFAngle.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFAngleList : public CMDFFList<CMDFFAngle>
{
public:
    CMDFFAngleList();
    virtual ~CMDFFAngleList();
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    void add(const CMDFFAngle& angle);
    void del(int index) { angles_.erase(angles_.begin() + index); }
    CMDFFAngle* get(int index) const { return angles_[index].get(); }
    int size() const { return (int)angles_.size(); }
    void empty() { angles_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string aAtom1, std::string atom2, std::string atom3) const;
    
    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFAngle* getRegisteredFFType(int index) const { return registeredForceFieldTypes_[index].get(); }

private:
    std::vector<std::shared_ptr<CMDFFAngle>> angles_;
    std::vector<std::shared_ptr<CMDFFAngle>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
