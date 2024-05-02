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

#include "MolecularSystemTools.h"

BEGIN_CUDA_COMPATIBLE()

void CMolecularSystemTools::removeDuplicateBondIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>* auxList) const
{
    int     atom1, atom2;
    int     numElements;
    bool    foundDuplicate;

    numElements = (int)atoms1.size();

    for(int i=0; i<numElements; i++)
    {
        atom1 = atoms1[i];
        atom2 = atoms2[i];

        for(int j=i+1; j<numElements; j++)
        {
            foundDuplicate = false;
            if(atom1 == atoms1[j])
            {
                if(atom2 == atoms2[j])
                {
                    foundDuplicate = true;
                }
            }

            if(atom1 == atoms2[j])
            {
                if(atom2 == atoms1[j])
                {
                    foundDuplicate = true;
                }
            }

            if(foundDuplicate)
            {
                atoms1.erase(atoms1.begin() + j);
                atoms2.erase(atoms2.begin() + j);
                if(auxList) auxList->erase(auxList->begin() + j);
                numElements--;
                j--;
            }
        }
    }
}

void CMolecularSystemTools::removeDuplicateAngleIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>* auxList) const
{
    int indexAtom1, indexAtom2, indexAtom3;
    int numElements;
    bool foundDuplicate;


    // Remove duplicate angles
    numElements = (int)atoms1.size();
    for(int i=0; i<numElements; i++)
    {
        indexAtom1 = atoms1[i];
        indexAtom2 = atoms2[i];
        indexAtom3 = atoms3[i];

        for(int j=i+1; j<numElements; j++)
        {
            foundDuplicate = false;
            if(indexAtom1 == atoms1[j])
            {
                if(indexAtom2 == atoms2[j])
                {
                    if(indexAtom3 == atoms3[j])
                    {
                        foundDuplicate = true;
                    }
                }
            }

            if(indexAtom1 == atoms3[j])
            {
                if(indexAtom2 == atoms2[j])
                {
                    if(indexAtom3 == atoms1[j])
                    {
                        foundDuplicate = true;
                    }
                }
            }

            if(foundDuplicate)
            {
                atoms1.erase(atoms1.begin() + j);
                atoms2.erase(atoms2.begin() + j);
                atoms3.erase(atoms3.begin() + j);
                if(auxList) auxList->erase(auxList->begin() + j);
                numElements--;
                j--;
            }
        }
    }


    // Remove angles with duplicate atoms
    numElements = (int)atoms1.size();
    for(int i=0; i<numElements; i++)
    {
        indexAtom1 = atoms1[i];
        indexAtom2 = atoms2[i];
        indexAtom3 = atoms3[i];

        foundDuplicate = false;
        if(indexAtom1 == indexAtom2) foundDuplicate = true;
        if(indexAtom1 == indexAtom3) foundDuplicate = true;
        if(indexAtom2 == indexAtom3) foundDuplicate = true;

        if(foundDuplicate)
        {
            atoms1.erase(atoms1.begin() + i);
            atoms2.erase(atoms2.begin() + i);
            atoms3.erase(atoms3.begin() + i);
            if(auxList) auxList->erase(auxList->begin() + i);
            numElements--;
            i--;
        }
    }
}

void CMolecularSystemTools::removeDuplicateDihIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4, std::vector<int>* auxList) const
{
    int indexAtom1, indexAtom2, indexAtom3, indexAtom4;
    int numElements;
    bool foundDuplicate;


    // Remove duplicate dihedrals
    numElements = (int)atoms1.size();
    for(int i=0; i<numElements; i++)
    {
        indexAtom1 = atoms1[i];
        indexAtom2 = atoms2[i];
        indexAtom3 = atoms3[i];
        indexAtom4 = atoms4[i];

        for(int j=i+1; j<numElements; j++)
        {
            foundDuplicate = false;
            if(indexAtom1 == atoms1[j])
            {
                if(indexAtom2 == atoms2[j])
                {
                    if(indexAtom3 == atoms3[j])
                    {
                        if(indexAtom4 == atoms4[j])
                        {
                            foundDuplicate = true;
                        }
                    }
                }
            }

            if(indexAtom1 == atoms4[j])
            {
                if(indexAtom2 == atoms3[j])
                {
                    if(indexAtom3 == atoms2[j])
                    {
                        if(indexAtom4 == atoms1[j])
                        {
                            foundDuplicate = true;
                        }
                    }
                }
            }

            if(foundDuplicate)
            {
                atoms1.erase(atoms1.begin() + j);
                atoms2.erase(atoms2.begin() + j);
                atoms3.erase(atoms3.begin() + j);
                atoms4.erase(atoms4.begin() + j);
                if(auxList) auxList->erase(auxList->begin() + j);
                numElements--;
                j--;
            }
        }
    }


    // Remove dihedrals with duplicate atoms
    numElements = (int)atoms1.size();
    for(int i=0; i<numElements; i++)
    {
        indexAtom1 = atoms1[i];
        indexAtom2 = atoms2[i];
        indexAtom3 = atoms3[i];
        indexAtom4 = atoms4[i];

        foundDuplicate = false;
        if(indexAtom1 == indexAtom2) foundDuplicate = true;
        if(indexAtom1 == indexAtom3) foundDuplicate = true;
        if(indexAtom1 == indexAtom4) foundDuplicate = true;
        if(indexAtom2 == indexAtom3) foundDuplicate = true;
        if(indexAtom2 == indexAtom4) foundDuplicate = true;
        if(indexAtom3 == indexAtom4) foundDuplicate = true;

        if(foundDuplicate)
        {
            atoms1.erase(atoms1.begin() + i);
            atoms2.erase(atoms2.begin() + i);
            atoms3.erase(atoms3.begin() + i);
            atoms4.erase(atoms4.begin() + i);
            if(auxList) auxList->erase(auxList->begin() + i);
            numElements--;
            i--;
        }
    }
}

void CMolecularSystemTools::genMolIndices(const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molIndices) const
{
    int nextMolIndex = 0;

    // Creates the MolIndices array. Each element will contain
    // a molecule ID (where the array index correspond to an atom index)
    molIndices.resize(bondDestIndices.size(), -1);

    for(int i=0; i<(int)bondDestIndices.size(); i++)
    {
        if(molIndices[i] != -1) continue;

        std::vector<int> aMolecule;
        getMoleculeConnectedToIndex(i, bondDestIndices, aMolecule);

        for(int j=0; j<(int)aMolecule.size(); j++)
            molIndices[aMolecule[j]] = nextMolIndex;

        nextMolIndex++;
    }
}

void CMolecularSystemTools::getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule) const
{
    std::map<int, int> visitationTracker;
    getMoleculeConnectedToIndex(index, bondDestIndices, molecule, visitationTracker);
}

void CMolecularSystemTools::getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule, std::map<int, int>& visitationTracker) const
{
    if(visitationTracker.count(index) == 1) return;

    molecule.emplace_back(index);
    visitationTracker[index] = 1;

    for(int i=0; i<(int)bondDestIndices[index].size(); i++)
    {
        int iNextIndex = bondDestIndices[index][i];
        if(iNextIndex >= (int)bondDestIndices.size()) continue;

        getMoleculeConnectedToIndex(iNextIndex, bondDestIndices, molecule, visitationTracker);
    }
}

END_CUDA_COMPATIBLE()
