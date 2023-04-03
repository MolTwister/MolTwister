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

#include "MolTwisterMDFFDihList.h"
#include "MolTwisterMDFFDih_Fourier4t.h"
#include "MolTwisterMDFFDih_Harm.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFDihList::CMDFFDihList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFDih_Fourier4t>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFDih_Harm>());
}

CMDFFDihList::~CMDFFDihList()
{
}

void CMDFFDihList::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFList<CMDFFDih>::serialize(io, saveToStream);

    // Note: we do not serialize registeredForceFieldTypes_, size this is always built in the
    // class constructor. Moreover, it is used to read from the stream
    if(saveToStream)
    {
        io << dihedrals_.size();
        for(std::shared_ptr<CMDFFDih> item : dihedrals_)
        {
            io << item->getFFType();
            item->serialize(io, saveToStream);
        }
    }
    else
    {
        size_t size;

        io >> size;
        dihedrals_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string ffType;
            io >> ffType;
            for(const std::shared_ptr<CMDFFDih>& item : registeredForceFieldTypes_)
            {
                if(item->getFFType() == ffType)
                {
                    dihedrals_[i] = item->createCopy();
                    dihedrals_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }
    }
}

void CMDFFDihList::add(const CMDFFDih& dih)
{    
    std::shared_ptr<std::vector<int>> indices = indexFromNames(dih.getAtomInBond(0), dih.getAtomInBond(1), dih.getAtomInBond(2), dih.getAtomInBond(3));
    appendToList(dihedrals_, *indices, dih);
}

std::shared_ptr<std::vector<int>> CMDFFDihList::indexFromNames(std::string atom1, std::string atom2, std::string atom3, std::string atom4) const
{
    auto indices = std::make_shared<std::vector<int>>();
    
    for(int i=0; i<(int)dihedrals_.size(); i++)
    {
        if(atom1 == dihedrals_[i]->getAtomInBond(0))
        {
            if(atom2 == dihedrals_[i]->getAtomInBond(1))
            {
                if(atom3 == dihedrals_[i]->getAtomInBond(2))
                {
                    if(atom4 == dihedrals_[i]->getAtomInBond(3))
                    {
                        indices->emplace_back(i);
                        continue;
                    }
                }
            }
        }
        
        if(atom1 == dihedrals_[i]->getAtomInBond(3))
        {
            if(atom2 == dihedrals_[i]->getAtomInBond(2))
            {
                if(atom3 == dihedrals_[i]->getAtomInBond(1))
                {
                    if(atom4 == dihedrals_[i]->getAtomInBond(0))
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
