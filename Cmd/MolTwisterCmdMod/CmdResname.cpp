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

#include "CmdResname.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdResname::getCmd()
{
    return "resname";
}

std::vector<std::string> CCmdResname::getCmdLineKeywords()
{
    return { "resname", "id", "var", "name", "resname", "molecule", "atomnames", "to" };
}

std::vector<std::string> CCmdResname::getCmdHelpLines()
{
    return {
                "resname id <atom index> to <resname>",
                "resname var <variable name> to <resname>",
                "resname name <atom ID> to <resname>",
                "resname resname <resname> to <resname>",
                "resname molecule <molecule index> to <resname>",
                "resname atomnames <atom ID1> <atom ID2> ... <atom IDn> to <resname>"
           };
}

std::string CCmdResname::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tModifies the resname of atoms given by\r\n";
    text+= "\t* the atom index, <atom index>\r\n";
    text+= "\t* the atom index contained in the variable <variable name>\r\n";
    text+= "\t* all atoms with atom ID = <atom ID> (e.g., O, H, C7)\r\n";
    text+= "\t* all atoms with resname = <resname>\r\n";
    text+= "\t* all atoms within molecular index = <molecule index>\r\n";
    text+= "\t* all atoms with atom IDs, <atom ID1> through <atom IDn>\r\n";
    text+= "\tto the resname given by <resname>.";

    return text;
}

std::string CCmdResname::execute(std::vector<std::string> arguments)
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
    else if(text1 == "resname")
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            text2 = state_->atoms_[i]->resname_;
            if(text2 == text1)
            {
                atomIndices.emplace_back(i);
                foundIndex = true;
            }
        }
    }
    else if(text1 == "molecule")
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);
        int molIndex = atoi(text1.data());

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->getMolIndex() == molIndex)
            {
                atomIndices.emplace_back(i);
                foundIndex = true;
            }
        }
    }
    else if(text1 == "atomnames")
    {
        int strLen;
        bool to;

        do
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);
            strLen = (int)text1.size();
            to = (text1 == "to") ? true : false;
            if((strLen > 0) && !to)
            {
                for(int i=0; i<state_->atoms_.size(); i++)
                {
                    std::string ID = state_->atoms_[i]->getID();
                    if(ID == text1)
                    {
                        atomIndices.emplace_back(i);
                        foundIndex = true;
                    }
                }
            }

        } while((strLen != 0) && !to);
        arg--;
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'id', 'var', 'name', 'resname', 'molecule' or 'atomnames'!";
        return lastError_;
    }

    if(foundIndex)
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 != "to")
        {
            lastError_= "Syntax Error: fifth argument should be 'to'!";
            return lastError_;
        }

        text1 = CASCIIUtility::getArg(arguments, arg++);

        for(int i=0; i<atomIndices.size(); i++)
        {
            int atomIndex = atomIndices[i];

            if(atomIndex < state_->atoms_.size())
            {
                atomPtr = state_->atoms_[atomIndex].get();
                if(atomPtr)
                {
                    atomPtr->resname_ = text1;
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
