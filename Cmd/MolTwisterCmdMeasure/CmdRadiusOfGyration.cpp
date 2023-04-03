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

#include "CmdRadiusOfGyration.h"
#include <functional>
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolTwisterStateTools.h"

namespace MolTwisterCmdMeasure
{
    std::string CCmdRadiusOfGyration::getCmd()
    {
        return "radiusofgyration";
    }

    std::vector<std::string> CCmdRadiusOfGyration::getCmdLineKeywords()
    {
        return { "radiusofgyration", "sel" };
    }

    std::vector<std::string> CCmdRadiusOfGyration::getCmdHelpLines()
    {
        return {
                    "radiusofgyration sel"
               };
    }

    std::string CCmdRadiusOfGyration::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the radius of gyration for the collection of atoms defined by\r\n";
        text+= "\t* the visual selection of atoms, by using the 'sel' keyword\r\n";
        text+= "\tIf 'sel' is not specified, 'sel' is assumed.\r\n";
        text+= "\r\n";
        text+= "\tNote that atomic masses must be loaded for this command to work!\r\n";
        text+= "\r\n";
        text+= "\tOutput:\r\n";
        text+= "\tRgyr = <radius of gyration>";

        return text;
    }

    std::string CCmdRadiusOfGyration::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text;
        double Rgyr = 0.0;
        int typeOfSelection = 0;

        // Type of selections possible are
        // 0: of selected atoms (default), 1: of selected atoms
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text.size() == 0)    typeOfSelection = 0;
        else if(text == "sel")  typeOfSelection = 0;
        else
        {
            lastError_ = "Error, type of selection not recognized!";
            return lastError_;
        }

        std::function<double(CMolTwisterState*, FILE*)> rgyrFromSel = [](CMolTwisterState* state, FILE* stdOut)
        {
            std::vector<int> atomIndices;

            for(int i=0; i<state->atoms_.size(); i++)
            {
                if(state->atoms_[i]->isSelected()) atomIndices.emplace_back(i);
            }

            C3DVector Rc = CMolTwisterStateTools(state, stdOut).getCenterOfMass(atomIndices, state->currentFrame_);

            double sum = 0.0;
            double sumM = 0.0;
            for(int i=0; i<state->atoms_.size(); i++)
            {
                if(state->atoms_[i]->isSelected())
                {
                    double m = state->atoms_[i]->m_;
                    C3DVector R = state->atoms_[i]->r_[state->currentFrame_];
                    C3DVector RmRc = R - Rc;
                    sum+= m * (RmRc * RmRc);
                    sumM+= m;
                }
            }

            return (sumM != 0.0) ? sqrt(sum / sumM) : -1.0;
        };

        if(typeOfSelection == 0)
        {
            Rgyr = rgyrFromSel(state_, stdOut_);
        }
        else
        {
            Rgyr = rgyrFromSel(state_, stdOut_);
        }

        fprintf(stdOut_, "\r\n\tRgyr = %.8f\r\n", Rgyr);

        return lastError_;
    }
}
