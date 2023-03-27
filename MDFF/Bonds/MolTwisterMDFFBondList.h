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
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFBond.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFBondList : public CMDFFList<CMDFFBond>
{
public:
    CMDFFBondList();
    virtual ~CMDFFBondList();
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    void add(const CMDFFBond& bond);
    void del(int index) { bonds_.erase(bonds_.begin() + index); }
    CMDFFBond* get(int index) const { return bonds_[index].get(); }
    int size() const { return (int)bonds_.size(); }
    void empty() { bonds_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2) const;

    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFBond* getRegisteredFFType(int index) const { return registeredForceFieldTypes_[index].get(); }
    
private:
    std::vector<std::shared_ptr<CMDFFBond>> bonds_;
    std::vector<std::shared_ptr<CMDFFBond>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
