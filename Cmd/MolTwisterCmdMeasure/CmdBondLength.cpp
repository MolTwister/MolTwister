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

#include "CmdBondLength.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"

namespace MolTwisterCmdMeasure
{
    std::string CCmdBondLength::getCmd()
    {
        return "bondlength";
    }

    std::vector<std::string> CCmdBondLength::getCmdLineKeywords()
    {
        return { "bondlength", "id", "var" };
    }

    std::vector<std::string> CCmdBondLength::getCmdHelpLines()
    {
        return {
                    "bondlength id <atom index 1> <atom index 2>",
                    "bondlength var <variable name>"
               };
    }

    std::string CCmdBondLength::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the length between two atoms by\r\n";
        text+= "\t* defining two atom indices directly by using the 'id' keyword\r\n";
        text+= "\t* obtaining atom indices from a variable by using the 'var' keyword\r\n";
        text+= "\r\n";
        text+= "\tOutput:\r\n";
        text+= "\tR(<atom index 1>, <atom index 2>) = <distance>";

        return text;
    }

    std::string CCmdBondLength::execute(std::vector<std::string> arguments)
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
                if((state_->getCurrFrameIndex() < atom1Ptr->r_.size()) && (state_->getCurrFrameIndex() < atom2Ptr->r_.size()))
                {
                    if(atom1Ptr && atom2Ptr)
                    {
                        double dDist = CMolecularTools(state_, stdOut_).measureDist(atom1Ptr, atom2Ptr, state_->getCurrFrameIndex());
                        fprintf(stdOut_, "\r\n\tR(%i, %i) = %g\r\n", atomIndex1, atomIndex2, dDist);
                    }
                    else
                    {
                        lastError_ = "Could not find requested atoms!";
                        return lastError_;
                    }
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
}
