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

#include "CmdAngle.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"

namespace MolTwisterCmdMod
{
    std::string CCmdAngle::getCmd()
    {
        return "angle";
    }

    std::vector<std::string> CCmdAngle::getCmdLineKeywords()
    {
        return { "angle", "id", "var" };
    }

    std::vector<std::string> CCmdAngle::getCmdHelpLines()
    {
        return {
                    "angle id <atom index 1> <atom index 2> <atom index 3> to <angle>",
                    "angle var <variable name> to <angle>"
               };
    }

    std::string CCmdAngle::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tModifies the angle, given by\r\n";
        text+= "\t* the three atom indices <atom index 1> to <atom index 3>\r\n";
        text+= "\t* the atom indices contained in the variable <variable name>\r\n";
        text+= "\tto the angle given by <angle>, which is given in degrees.";

        return text;
    }

    std::string CCmdAngle::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text1, text2;
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        CAtom* atom3Ptr = nullptr;
        bool foundIndex = false;
        int atomIndex1=0, atomIndex2=0, atomIndex3=0;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 == "id")
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndex1 = atoi(text1.data());
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndex2 = atoi(text1.data());
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndex3 = atoi(text1.data());
            foundIndex = true;
        }
        else if(text1 == "var")
        {
            int varIndex;

            text1 = CASCIIUtility::getArg(arguments, arg++);
            CVar* varPtr = state_->getVariable(text1.data(), varIndex);
            if(varPtr && (varPtr->getType() == CVar::typeAngle))
            {
                CVarAngle* p = (CVarAngle*)varPtr;

                atomIndex1 = p->atomIndex1_;
                atomIndex2 = p->atomIndex2_;
                atomIndex3 = p->atomIndex3_;
                foundIndex = true;
            }
            else
            {
                lastError_ = "Variable missing or not of type 'angle'!";
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
            if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size())
               && (atomIndex3 < state_->atoms_.size()))
            {
                atom1Ptr = state_->atoms_[atomIndex1].get();
                atom2Ptr = state_->atoms_[atomIndex2].get();
                atom3Ptr = state_->atoms_[atomIndex3].get();
                if(atom1Ptr && atom2Ptr && atom3Ptr)
                {
                    double angle;

                    text1 = CASCIIUtility::getArg(arguments, arg++);
                    text2 = text1;

                    text1 = CASCIIUtility::getArg(arguments, arg++);
                    angle = atof(text1.data()) * M_PI / 180.0;

                    if((state_->currentFrame_ < atom1Ptr->r_.size()) && (state_->currentFrame_ < atom2Ptr->r_.size())
                       && (state_->currentFrame_ < atom3Ptr->r_.size()))
                    {
                        if(text2 == "to")
                        {
                            CMolecularTools::modAngleTo(atom1Ptr, atom2Ptr, atom3Ptr, angle, state_->currentFrame_);
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
