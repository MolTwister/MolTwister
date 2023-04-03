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

#include "MolTwisterMDFFAngleList.h"
#include "MolTwisterMDFFAngle_Harm.h"
#include "MolTwisterMDFFAngle_Class2.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFAngleList::CMDFFAngleList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFAngle_Harm>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFAngle_Class2>());
}

CMDFFAngleList::~CMDFFAngleList()
{
}

void CMDFFAngleList::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFList<CMDFFAngle>::serialize(io, saveToStream);

    // Note: we do not serialize registeredForceFieldTypes_, size this is always built in the
    // class constructor. Moreover, it is used to read from the stream
    if(saveToStream)
    {
        io << angles_.size();
        for(std::shared_ptr<CMDFFAngle> item : angles_)
        {
            io << item->getFFType();
            item->serialize(io, saveToStream);
        }
    }
    else
    {
        size_t size;

        io >> size;
        angles_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string ffType;
            io >> ffType;
            for(const std::shared_ptr<CMDFFAngle>& item : registeredForceFieldTypes_)
            {
                if(item->getFFType() == ffType)
                {
                    angles_[i] = item->createCopy();
                    angles_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }
    }
}

void CMDFFAngleList::add(const CMDFFAngle& angle)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(angle.getAtomInBond(0), angle.getAtomInBond(1), angle.getAtomInBond(2));
    appendToList(angles_, *indices, angle);
}

std::shared_ptr<std::vector<int>> CMDFFAngleList::indexFromNames(std::string aAtom1, std::string atom2, std::string atom3) const
{
    auto indices = std::make_shared<std::vector<int>>();

    for(int i=0; i<(int)angles_.size(); i++)
    {
        if(aAtom1 == angles_[i]->getAtomInBond(0))
        {
            if(atom2 == angles_[i]->getAtomInBond(1))
            {
                if(atom3 == angles_[i]->getAtomInBond(2))
                {
                    indices->emplace_back(i);
                    continue;
                }
            }
        }
        
        if(angles_[i]->isAngleABCEqualToCBA())
        {
            if(aAtom1 == angles_[i]->getAtomInBond(2))
            {
                if(atom2 == angles_[i]->getAtomInBond(1))
                {
                    if(atom3 == angles_[i]->getAtomInBond(0))
                    {
                        indices->emplace_back(i);
                        continue;
                    }
                }
            }
        }
    }

    return indices;
}

END_CUDA_COMPATIBLE()
