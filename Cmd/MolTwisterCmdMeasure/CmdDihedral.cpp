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

#include "CmdDihedral.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"
#include "../Tools/ParsingTools.h"

namespace MolTwisterCmdMeasure
{
    std::string CCmdDihedral::getCmd()
    {
        return "dihedral";
    }

    std::vector<std::string> CCmdDihedral::getCmdLineKeywords()
    {
        return { "dihedral", "id", "var" };
    }

    std::vector<std::string> CCmdDihedral::getCmdHelpLines()
    {
        return {
                    "dihedral <dihedral>"
               };
    }

    std::string CCmdDihedral::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the dihedral angle, given by the <dihedral> parameter. This\r\n";
        text+= "\tparameter specifies the dihedral to rotate and can be formatted\r\n";
        text+= "\tas shown below.\r\n";
        text+= "\t* id <atom index 1> <atom index 2> <atom index 3> <atom index 4>\r\n";
        text+= "\t* var <variable name>\r\n";
        text+= "\r\n";
        text+= "\tOutput:\r\n";
        text+= "\tPhi(<dih. index 1>, <dih. index 2>, <dih. index 3>, <dih. index 4>) = <angle>\r\n";
        text+= "\twhere <angle> is given in degrees.";

        return text;
    }

    std::string CCmdDihedral::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        CAtom* atom3Ptr = nullptr;
        CAtom* atom4Ptr = nullptr;
        int atomIndices[4];

        std::pair<bool, std::string> retVal = CParsingTools(state_, stdOut_).retrieveDihedralAtoms(arguments, arg, &atom1Ptr, &atom2Ptr, &atom3Ptr, &atom4Ptr, atomIndices);
        if(retVal.first)
        {
            if((state_->getCurrFrameIndex() < atom1Ptr->r_.size()) && (state_->getCurrFrameIndex() < atom2Ptr->r_.size())
               && (state_->getCurrFrameIndex() < atom3Ptr->r_.size()) && (state_->getCurrFrameIndex() < atom4Ptr->r_.size()))
            {
                if(atom1Ptr && atom2Ptr && atom3Ptr && atom4Ptr)
                {
                    double dPhi = CMolecularTools(state_, stdOut_).measureDihedral(atom1Ptr, atom2Ptr, atom3Ptr, atom4Ptr, state_->getCurrFrameIndex());
                    fprintf(stdOut_, "\r\nPhi(%i, %i, %i, %i) = %g\r\n", atomIndices[0], atomIndices[1], atomIndices[2], atomIndices[3], dPhi * 180.0/M_PI);
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
            lastError_ = retVal.second;
        }

        return lastError_;
    }
}
