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

#include "MolTwisterMDFFBondList.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFBondList::CMDFFBondList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_Harm>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_Morse>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_LJC>());
}

CMDFFBondList::~CMDFFBondList()
{
}

void CMDFFBondList::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFList<CMDFFBond>::serialize(io, saveToStream);

    // Note: we do not serialize registeredForceFieldTypes_, size this is always built in the
    // class constructor. Moreover, it is used to read from the stream
    if(saveToStream)
    {
        io << bonds_.size();
        for(std::shared_ptr<CMDFFBond> item : bonds_)
        {
            io << item->getFFType();
            item->serialize(io, saveToStream);
        }
    }
    else
    {
        size_t size;

        io >> size;
        bonds_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string ffType;
            io >> ffType;
            for(std::shared_ptr<CMDFFBond> item : registeredForceFieldTypes_)
            {
                if(item->getFFType() == ffType)
                {
                    bonds_[i] = item->createCopy();
                    bonds_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }
    }
}

void CMDFFBondList::add(const CMDFFBond& bond)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(bond.getAtomInBond(0), bond.getAtomInBond(1));
    appendToList(bonds_, *indices, bond);
}

std::shared_ptr<std::vector<int>> CMDFFBondList::indexFromNames(std::string atom1, std::string atom2) const
{
    auto indices = std::make_shared<std::vector<int>>();
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(atom1 == bonds_[i]->getAtomInBond(0))
        {
            if(atom2 == bonds_[i]->getAtomInBond(1))
            {
                indices->emplace_back(i);
                continue;
            }
        }
        
        if(atom1 == bonds_[i]->getAtomInBond(1))
        {
            if(atom2 == bonds_[i]->getAtomInBond(0))
            {
                indices->emplace_back(i);
                continue;
            }
        }
    }

    return indices;
}

END_CUDA_COMPATIBLE()
