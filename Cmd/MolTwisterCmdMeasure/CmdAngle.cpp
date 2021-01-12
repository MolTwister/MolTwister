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

namespace MolTwisterCmdMeasure
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
                    "angle id <atom index 1> <atom index 2> <atom index 3>",
                    "angle var <variable name>"
               };
    }

    std::string CCmdAngle::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the angle between three atoms (in degrees) by\r\n";
        text+= "\t* defining three atom indices directly by using the 'id' keyword\r\n";
        text+= "\t* obtaining atom indices from a variable by using the 'var' keyword\r\n";
        text+= "\r\n";
        text+= "\tOutput:\r\n";
        text+= "\tTheta(<atom index 1>, <atom index 2>, <atom index 3>) = <angle in degrees>";

        return text;
    }

    std::string CCmdAngle::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text;
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        CAtom* atom3Ptr = nullptr;
        bool foundIndex = false;
        int atomIndex1=0, atomIndex2=0, atomIndex3=0;


        text = CASCIIUtility::getArg(arguments, arg++);
        if(text == "id")
        {
            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex1 = atoi(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex2 = atoi(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex3 = atoi(text.data());
            foundIndex = true;
        }
        else if(text == "var")
        {
            int varIndex;

            text = CASCIIUtility::getArg(arguments, arg++);
            CVar* varPtr = state_->getVariable(text.data(), varIndex);
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
                if((state_->getCurrFrameIndex() < atom1Ptr->r_.size()) && (state_->getCurrFrameIndex() < atom2Ptr->r_.size())
                   && (state_->getCurrFrameIndex() < atom3Ptr->r_.size()))
                {
                    if(atom1Ptr && atom2Ptr && atom3Ptr)
                    {
                        double theta = CMolecularTools(state_, stdOut_).measureAngle(atom1Ptr, atom2Ptr, atom3Ptr, state_->getCurrFrameIndex());
                        fprintf(stdOut_, "\r\nTheta(%i, %i, %i) = %g\r\n", atomIndex1, atomIndex2, atomIndex3, theta * 180.0/M_PI);
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
