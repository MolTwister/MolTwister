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

#include "CmdBondSep.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdBondSep::getCmd()
{
    return "bondsep";
}

std::vector<std::string> CCmdBondSep::getCmdLineKeywords()
{
    return { "bondsep", "id", "var" };
}

std::vector<std::string> CCmdBondSep::getCmdHelpLines()
{
    return {
                "bondsep id <atom index 1> <atom index 2>",
                "bondsep var <variable name>"
           };
}

std::string CCmdBondSep::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures the number of bondsbetween two atoms by\r\n";
    text+= "\t* defining two atom indices directly by using the 'id' keyword\r\n";
    text+= "\t* obtaining atom indices from a variable by using the 'var' keyword\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t* BondSep(<atom index 1>, <atom index 2>) = <distance>\r\n";
    text+= "\tor\r\n";
    text+= "\t* BondSep(<atom index 1>, <atom index 2>) > 4\r\n";
    text+= "\tif the number of bonds exceeds 4.";

    return text;
}

std::string CCmdBondSep::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CAtom* atom1Ptr = nullptr;
    CAtom* atom2Ptr = nullptr;
    bool foundIndex = false;
    int atomIndex1=0, atomIndex2=0;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndex1 = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndex2 = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;

        text = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeBond))
        {
            CVarBond* p = (CVarBond*)varPtr;

            atomIndex1 = p->atomIndex1_;
            atomIndex2 = p->atomIndex2_;
            foundIndex = true;
        }
        else
        {
            lastError_ = "Variable missing or not of type 'bond'!";
            return lastError_;
        }
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'id' or 'var'!";
        return lastError_;
    }

    if(foundIndex)
    {
        if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size()))
        {
            atom1Ptr = state_->atoms_[atomIndex1].get();
            atom2Ptr = state_->atoms_[atomIndex2].get();
            if(atom1Ptr && atom2Ptr)
            {
                int iDist = atom1Ptr->getBondSepTo(atom2Ptr);
                if(iDist == -1) fprintf(stdOut_, "\r\nBondSep(%i, %i) > 4\r\n", atomIndex1, atomIndex2);
                else            fprintf(stdOut_, "\r\nBondSep(%i, %i) = %i\r\n", atomIndex1, atomIndex2, iDist);

            }
            else
            {
                lastError_ = "Could not find requested atoms!";
                return lastError_;
            }
        }
        else
        {
            lastError_ = "Invalid atom indices!";
            return lastError_;
        }
    }

    return lastError_;
}
