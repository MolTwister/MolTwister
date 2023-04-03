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
#include <string.h>
#include <vector>
#include <memory>
#include "Utilities/Serializer.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

template<class ArrayType>
class CMDFFList
{
public:
    CMDFFList() {}
    
public:
    virtual void serialize(CSerializer&, bool saveToStream)
    {
        if(saveToStream)
        {
            // Plaeholder for future serialization!
        }
        else
        {
            // Plaeholder for future serialization!
        }
    }

    void appendToList(std::vector<std::shared_ptr<ArrayType>>& list, const std::vector<int>& indicesOfSimilarEntries, const ArrayType& itemToAppend) const
    {
        if(indicesOfSimilarEntries.size() == 0)
            list.emplace_back(itemToAppend.createCopy());
        else
        {
            int indexOfExistingItem = -1;
            for(int i=0; i<indicesOfSimilarEntries.size(); i++)
            {
                if(itemToAppend.getFFType() == list[indicesOfSimilarEntries[i]]->getFFType())
                {
                    indexOfExistingItem = indicesOfSimilarEntries[i];
                    break;
                }
            }
            
            if(indexOfExistingItem != -1)
            {
                list.erase(list.begin() + indexOfExistingItem);
                list.insert(list.begin() + indexOfExistingItem, itemToAppend.createCopy());
            }
            else
            {
                list.emplace_back(itemToAppend.createCopy());
            }
        }
    }
};

END_CUDA_COMPATIBLE()
