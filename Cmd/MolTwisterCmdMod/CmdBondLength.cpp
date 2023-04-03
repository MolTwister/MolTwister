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

namespace MolTwisterCmdMod
{
    std::string CCmdBondLength::getCmd()
    {
        return "bondlength";
    }

    std::vector<std::string> CCmdBondLength::getCmdLineKeywords()
    {
        return { "bondlength", "id", "var", "to" };
    }

    std::vector<std::string> CCmdBondLength::getCmdHelpLines()
    {
        return {
                    "bondlength id <atom index 1> <atom index 2> to <dist>",
                    "bondlength var <variable name> to <dist>"
               };
    }

    std::string CCmdBondLength::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tModifies the distance between two atoms, typically a bond, given by\r\n";
        text+= "\t* the two atom indices <atom index 1> and <atom index 3>\r\n";
        text+= "\t* the atom indices contained in the variable <variable name>\r\n";
        text+= "\tto the distance given by <dist>.";

        return text;
    }

    std::string CCmdBondLength::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text1, text2;
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        bool foundIndex = false;
        int atomIndex1=0, atomIndex2=0;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 == "id")
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndex1 = atoi(text1.data());
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndex2 = atoi(text1.data());
            foundIndex = true;
        }
        else if(text1 == "var")
        {
            int varIndex;

            text1 = CASCIIUtility::getArg(arguments, arg++);
            CVar* varPtr = state_->getVariable(text1.data(), varIndex);
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
                    double dist;

                    text1 = CASCIIUtility::getArg(arguments, arg++);
                    text2 = text1;

                    text1 = CASCIIUtility::getArg(arguments, arg++);
                    dist = atof(text1.data());

                    if((state_->currentFrame_ < atom1Ptr->r_.size()) && (state_->currentFrame_ < atom2Ptr->r_.size()))
                    {
                        if(text2 == "to")
                        {
                            CMolecularTools::modBondLengthTo(atom1Ptr, atom2Ptr, dist, state_->currentFrame_);
                        }
                        else
                        {
                            lastError_ = "Syntax Error: fifth argument should be 'to'!";
                            return lastError_;
                        }
                    }
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
}
