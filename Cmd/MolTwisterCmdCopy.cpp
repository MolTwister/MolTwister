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

#include <iostream>
#include <vector>
#include "MolTwisterCmdCopy.h"

void CCmdCopy::onAddKeywords()
{
    addKeyword("copy");
    addKeyword("atoms");
    addKeyword("names");
    addKeyword("resnames");
    addKeyword("all");
    addKeyword("sel");
}

std::string CCmdCopy::getHelpString() const
{
    std::string text;
    
    text+= "\tUsage: copy <specifier>\r\n";
    text+= "\r\n";
    text+= "\tMake several copies of atoms, or collection of atoms. The allowed specifiers are:\r\n";
    text+= "\r\n";
    text+= "\t        * atoms <nx> <ny> <nz> <dx> <dy> <dz> <index1> ... <indexN>        : Copy <ni> times\r\n";
    text+= "\t          with distance <di> between each copy, where i is x,y or z.\r\n";
    text+= "\t        * names <nx> <ny> <nz> <dx> <dy> <dz> <name1> ... <nameN>          : Copy based on names.\r\n";
    text+= "\t        * all <nx> <ny> <nz> <dx> <dy> <dz>                                : Copy all.\r\n";
    text+= "\t        * sel <nx> <ny> <nz> <dx> <dy> <dz>                                : Copy only selected.\r\n";
    text+= "\t        * resnames <nx> <ny> <nz> <dx> <dy> <dz> <resname1> ... <resnameN> : Copy based on resnames.";
    
    return text;
}

void CCmdCopy::execute(std::string commandLine)
{
    std::vector<std::string> objectsToCopy;
    std::vector<CAtom*> atomsToCopy;
    std::string text, commandString;
    double dX, dY, dZ;
    int arg = 1, length, count = 0;
    int nX, nY, nZ;
    double sX, sY, sZ;
    
    if(!state_) return;
    
    commandString = CASCIIUtility::getWord(commandLine, arg++);
    if((commandString != "atoms") &&
       (commandString != "names") &&
       (commandString != "all") &&
       (commandString != "sel") &&
       (commandString != "resnames"))
    {
        printf("Syntax Error: First argument should be the type of object to copy!");
        return;
    }
    
    // Get amoount of copies to do and distances
    text = CASCIIUtility::getWord(commandLine, arg++);
    nX = atoi(text.data());
    text = CASCIIUtility::getWord(commandLine, arg++);
    nY = atoi(text.data());
    text = CASCIIUtility::getWord(commandLine, arg++);
    nZ = atoi(text.data());

    text = CASCIIUtility::getWord(commandLine, arg++);
    dX = atof(text.data());
    text = CASCIIUtility::getWord(commandLine, arg++);
    dY = atof(text.data());
    text = CASCIIUtility::getWord(commandLine, arg++);
    dZ = atof(text.data());
    
    // Get list of objects to copy
    do
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        length = (int)text.size();
        count++;
        
        if(length)
        {
            objectsToCopy.emplace_back(text);
        }
    
    } while(length && (count < 200));

    // Convert objects to proper atom pointers and produce a list: aAtoms to copy,
    if(commandString == "atoms")
    {
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            int iIndex = atoi(objectsToCopy[i].data());
            if((iIndex < 0) || (iIndex >= state_->atoms_.size()))
            {
                printf("Error: found no atom with index %i, aborted copy!", iIndex);
                return;
            }

            CAtom* atom = state_->atoms_[iIndex].get();
            if(!atom)
            {
                printf("Error: found empty atom with index %i, aborted copy!", iIndex);
                return;
            }
        }
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            int iIndex = atoi(objectsToCopy[i].data());
            atomsToCopy.emplace_back(state_->atoms_[iIndex].get());
        }
    }
    else if(commandString == "names")
    {
        std::vector<int> atoms;
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            state_->getAtomsWithID(objectsToCopy[i], atoms);

            for(int j=0; j<atoms.size(); j++)
            {
                atomsToCopy.emplace_back(state_->atoms_[atoms[j]].get());
            }
        }
    }
    else if(commandString == "all")
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            atomsToCopy.emplace_back(state_->atoms_[j].get());
        }
    }
    else if(commandString == "sel")
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            if(state_->atoms_[j]->isSelected())
                atomsToCopy.emplace_back(state_->atoms_[j].get());
        }
    }
    else if(commandString == "resnames")
    {
        std::vector<int> atoms;
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            state_->getAtomsWithResname(objectsToCopy[i], atoms);
            
            for(int j=0; j<atoms.size(); j++)
            {
                atomsToCopy.emplace_back(state_->atoms_[atoms[j]].get());
            }
        }
    }
    
    // Copy the atoms, keep only intrinsic atom properties in copy
    for(int i=0; i<atomsToCopy.size(); i++)
    {
        CAtom* atomPtr = atomsToCopy[i];
        
        sX = atomPtr->r_[0].x_;
        sY = atomPtr->r_[0].y_;
        sZ = atomPtr->r_[0].z_;
        
        CAtom atom;
        atom.copyIntrinsicAtomProperties(*atomPtr);
        
        for(int iX=0; iX<nX; iX++)
        {
            for(int iY=0; iY<nY; iY++)
            {
                for(int iZ=0; iZ<nZ; iZ++)
                {
                    // Make sure we do not make a copy of the atom that is already there
                    if((iX == 0) && (iY == 0) && (iZ == 0)) continue;
                    
                    // Copy to the new position away from pos. of existing atom
                    atom.r_[0].x_ = sX + double(iX)*dX;
                    atom.r_[0].y_ = sY + double(iY)*dY;
                    atom.r_[0].z_ = sZ + double(iZ)*dZ;
                    state_->addAtom(atom);
                }
            }
        }
    }

    // Update 3D view
    if(state_->view3D_)
    {
        state_->view3D_->requestUpdate(false);
    }
}
