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

#include "CmdMass.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdMass::getCmd()
{
    return "mass";
}

std::vector<std::string> CCmdMass::getCmdLineKeywords()
{
    return { "mass", "id", "var", "name" "to" };
}

std::vector<std::string> CCmdMass::getCmdHelpLines()
{
    return {
                "mass id <atom index> to <atomic mass>",
                "mass var <variable name> to <atomic mass>",
                "mass name <atom ID> to <atomic mass>"
           };
}

std::string CCmdMass::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tModifies the atomic mass of atoms given by\r\n";
    text+= "\t* the atom index, <atom index>\r\n";
    text+= "\t* the atom index contained in the variable <variable name>\r\n";
    text+= "\t* all atoms with atom ID = <atom ID> (e.g., O, H, C7)\r\n";
    text+= "\tto the mass given by <atomic mass>.";

    return text;
}

std::string CCmdMass::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text1, text2;
    CAtom* atomPtr = nullptr;
    bool foundIndex = false;
    std::vector<int> atomIndices;

    text1 = CASCIIUtility::getArg(arguments, arg++);
    if(text1 == "id")
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);
        atomIndices.emplace_back(atoi(text1.data()));
        foundIndex = true;
    }
    else if(text1 == "var")
    {
        int varIndex;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text1.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeAtom))
        {
            CVarAtom* p = (CVarAtom*)varPtr;

            atomIndices.emplace_back(p->atomIndex_);
            foundIndex = true;
        }
        else
        {
            lastError_ = "Variable missing or not of type 'atom'!";
            return lastError_;
        }
    }
    else if(text1 == "name")
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            text2 = state_->atoms_[i]->getID();
            if(text2 == text1)
            {
                atomIndices.emplace_back(i);
                foundIndex = true;
            }
        }
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'id', 'var' or 'name'!";
        return lastError_;
    }

    if(foundIndex)
    {
        double m;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 != "to")
        {
            lastError_ = "Syntax Error: fifth argument should be 'to'!";
            return lastError_;
        }

        text1 = CASCIIUtility::getArg(arguments, arg++);
        m = atof(text1.data());

        for(int i=0; i<atomIndices.size(); i++)
        {
            int atomIndex = atomIndices[i];

            if(atomIndex < state_->atoms_.size())
            {
                atomPtr = state_->atoms_[atomIndex].get();
                if(atomPtr)
                {
                    atomPtr->m_ = m;
                }
                else
                {
                    lastError_ = "Could not find requested atom!";
                    return lastError_;
                }
            }
            else
            {
                lastError_ = "Invalid atom index!";
                return lastError_;
            }
        }
    }

    return lastError_;
}
