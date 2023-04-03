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

#include "CmdCoulombPotential.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ParsingTools.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/MolecularTools.h"

std::string CCmdCoulombPotential::getCmd()
{
    return "coulombpotential";
}

std::vector<std::string> CCmdCoulombPotential::getCmdLineKeywords()
{
    return { "coulombpotential", "single", "dihedralrot", "rot", "id", "var", "at" };
}

std::vector<std::string> CCmdCoulombPotential::getCmdHelpLines()
{
    return {
                "coulombpotential single at <x> <y> <z>",
                "coulombpotential dihedralrot at <x> <y> <z> <dihedral> rot <angle start> <angle end> <angular step size>"
           };
}

std::string CCmdCoulombPotential::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures the Coulomb potential, at a given point, of a test charge, +1,\r\n";
    text+= "\tusing one of the following choices\r\n";
    text+= "\t* 'single':      find coulomb energy at <x> <y> <z>, with output 'Utot =\r\n";
    text+= "\t                  <energy> kJ/C'.\r\n";
    text+= "\t* 'dihedralrot': find coulomb energy at <x> <y> <z> for several dihedral\r\n";
    text+= "\t                 rotations between <angle start> and <angle end> with steps\r\n";
    text+= "\t                 <angular step size> (in degrees), where the rotation is\r\n";
    text+= "\t                 applied to the specified dihedral, <dihedral>.\r\n";
    text+= "\r\n";
    text+= "\tThe <dihedral> parameter specifies the dihedral to rotate and can be formatted\r\n";
    text+= "\tas shown below.\r\n";
    text+= "\t* id <atom index 1> <atom index 2> <atom index 3> <atom index 4>\r\n";
    text+= "\t* var <variable name>";

    return text;
}

std::string CCmdCoulombPotential::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CMolTwisterStateTools stateTools(state_, stdOut_);
    std::string text;
    int path = -1;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "single")
    {
        path = 1;
    }
    else if(text == "dihedralrot")
    {
        path = 2;
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should indicate for example single point, dihedral rotation, etc.!";
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "at")
    {
        C3DVector at;

        text = CASCIIUtility::getArg(arguments, arg++);
        at.x_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        at.y_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        at.z_ = atof(text.data());

        if(path == 1)
        {
            fprintf(stdOut_, "\r\nUtot = %g kJ/C\r\n", stateTools.measureCoulombPotential(at));
        }
        else if(path == 2)
        {
            CAtom* atom1Ptr = nullptr;
            CAtom* atom2Ptr = nullptr;
            CAtom* atom3Ptr = nullptr;
            CAtom* atom4Ptr = nullptr;
            int atomIndices[4];

            std::pair<bool, std::string> retVal = CParsingTools(state_, stdOut_).retrieveDihedralAtoms(arguments, arg, &atom1Ptr, &atom2Ptr, &atom3Ptr, &atom4Ptr, atomIndices);
            if(retVal.first)
            {
                double startAngle, endAngle, step;

                text = CASCIIUtility::getArg(arguments, arg++);
                if(text == "rot")
                {
                    text = CASCIIUtility::getArg(arguments, arg++);
                    startAngle = atof(text.data());
                    text = CASCIIUtility::getArg(arguments, arg++);
                    endAngle = atof(text.data());
                    text = CASCIIUtility::getArg(arguments, arg++);
                    step = atof(text.data());

                    fprintf(stdOut_, "\r\n\t%-15s%-15s\r\n\t------------------------------------\r\n", "Angle [Deg]", "Etot [kJ/mol]");
                    if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                    {
                        for(double angle=startAngle; angle<=endAngle; angle+=step)
                        {
                            CMolecularTools::modDihedralTo(atom1Ptr, atom2Ptr, atom3Ptr, atom4Ptr, angle * M_PI/180.0, state_->getCurrFrameIndex());
                            fprintf(stdOut_, "\t% -15g% -15g\r\n", angle, stateTools.measureCoulombPotential(at));
                        }
                        state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                    }
                }
                else
                {
                    lastError_ = "Syntax Error: Fifth argument should be 'rot'!";
                    return  lastError_;
                }
            }
            else
            {
                lastError_ = retVal.second;
                return lastError_;
            }
        }
        else
        {
            lastError_ = "Error: Unrecognized parameters!";
            return lastError_;
        }
    }
    else
    {
        lastError_ = "Syntax Error: Fourth argument should be 'at'!";
    }

    return lastError_;
}
