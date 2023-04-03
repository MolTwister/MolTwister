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

#include "MolTwisterMDFFNonBondedList.h"
#include "MolTwisterMDFFNonBonded_LJ.h"
#include "MolTwisterMDFFNonBonded_LJ1208.h"
#include "MolTwisterMDFFNonBonded_Buck.h"
#include "MolTwisterMDFFNonBonded_LJBuck.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFNonBondedList::CMDFFNonBondedList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJ>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJ1208>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_Buck>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJBuck>());
}

CMDFFNonBondedList::~CMDFFNonBondedList()
{
}

void CMDFFNonBondedList::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFList<CMDFFNonBonded>::serialize(io, saveToStream);

    // Note: we do not serialize registeredForceFieldTypes_, size this is always built in the
    // class constructor. Moreover, it is used to read from the stream
    if(saveToStream)
    {
        io << nonBonded_.size();
        for(std::shared_ptr<CMDFFNonBonded> item : nonBonded_)
        {
            io << item->getFFType();
            item->serialize(io, saveToStream);
        }
    }
    else
    {
        size_t size;

        io >> size;
        nonBonded_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string ffType;
            io >> ffType;
            for(const std::shared_ptr<CMDFFNonBonded>& item : registeredForceFieldTypes_)
            {
                if(item->getFFType() == ffType)
                {
                    nonBonded_[i] = item->createCopy();
                    nonBonded_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }
    }
}

void CMDFFNonBondedList::add(const CMDFFNonBonded& nonBonded)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(nonBonded.getAtomInBond(0), nonBonded.getAtomInBond(1));
    appendToList(nonBonded_, *indices, nonBonded);
}

std::shared_ptr<std::vector<int>> CMDFFNonBondedList::indexFromNames(std::string atom1, std::string atom2)
{
    auto indices = std::make_shared<std::vector<int>>();

    for(int i=0; i<(int)nonBonded_.size(); i++)
    {
        if(atom1 == nonBonded_[i]->getAtomInBond(0))
        {
            if(atom2 == nonBonded_[i]->getAtomInBond(1))
            {
                indices->emplace_back(i);
                continue;
            }
        }

        if(atom1 == nonBonded_[i]->getAtomInBond(1))
        {
            if(atom2 == nonBonded_[i]->getAtomInBond(0))
            {
                indices->emplace_back(i);
                continue;
            }
        }
    }

    return indices;
}

END_CUDA_COMPATIBLE()
