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

#include "CmdCenter.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdCenter::getCmd()
{
    return "center";
}

std::vector<std::string> CCmdCenter::getCmdLineKeywords()
{
    return { "center", "sel" };
}

std::vector<std::string> CCmdCenter::getCmdHelpLines()
{
    return {
                "center sel"
           };
}

std::string CCmdCenter::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures the geometric center of the collection of atoms defined by\r\n";
    text+= "\t* the visual selection of atoms, by using the 'sel' keyword\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t(<x_center>, <y_center>, <z_center>)";

    return text;
}

std::string CCmdCenter::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;

    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "sel")
    {
        C3DVector v;
        int N = 0;
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->isSelected())
            {
                if(!(state_->currentFrame_ < state_->atoms_[i]->r_.size())) continue;
                v+= state_->atoms_[i]->r_[state_->currentFrame_];
                N++;
            }
        }

        if(N > 0)
        {
            fprintf(stdOut_, "\r\n\t(%.4f, %.4f, %.4f)\r\n", v.x_ / double(N), v.y_ / double(N), v.z_ / double(N));
        }

        return lastError_;
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'sel'!";
    }

    return lastError_;
}
