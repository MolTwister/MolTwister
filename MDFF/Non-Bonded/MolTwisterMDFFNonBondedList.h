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
#include "MolTwisterMDFFNonBonded.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFNonBondedList : public CMDFFList<CMDFFNonBonded>
{
public:
    CMDFFNonBondedList();
    virtual ~CMDFFNonBondedList();
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    void add(const CMDFFNonBonded& nonBonded);
    void del(int index) { nonBonded_.erase(nonBonded_.begin() + index); }
    CMDFFNonBonded* get(int index) const { return nonBonded_[index].get(); }
    int size() const { return (int)nonBonded_.size(); }
    void empty() { nonBonded_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2);
    
    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFNonBonded* getRegisteredFFType(int iIndex) const { return registeredForceFieldTypes_[iIndex].get(); }

private:
    std::vector<std::shared_ptr<CMDFFNonBonded>> nonBonded_;
    std::vector<std::shared_ptr<CMDFFNonBonded>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
