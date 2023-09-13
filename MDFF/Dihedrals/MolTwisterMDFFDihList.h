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
#include "MolTwisterMDFFDih.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFDihList : public CMDFFList<CMDFFDih>
{
public:
    CMDFFDihList();
    virtual ~CMDFFDihList();
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    void add(const CMDFFDih& dih);
    void del(int index) { dihedrals_.erase(dihedrals_.begin() + index); }
    CMDFFDih* get(int index) const { return dihedrals_[index].get(); }
    int size() const { return (int)dihedrals_.size(); }
    void empty() { dihedrals_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2, std::string atom3, std::string atom4) const;

    int getNumRegisteredFFTypes() { return (int)registeredForceFieldTypes_.size(); }
    CMDFFDih* getRegisteredFFType(int index) { return registeredForceFieldTypes_[index].get(); }
    
private:
    std::vector<std::shared_ptr<CMDFFDih>> dihedrals_;
    std::vector<std::shared_ptr<CMDFFDih>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
