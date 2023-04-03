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

#include "CmdCount.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdCount::getCmd()
{
    return "count";
}

std::vector<std::string> CCmdCount::getCmdLineKeywords()
{
    return { "count", "sel", "tot" };
}

std::vector<std::string> CCmdCount::getCmdHelpLines()
{
    return {
                "count sel",
                "count tot"
           };
}

std::string CCmdCount::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCounts number of atoms. Either the total count can be measured, using 'tot',\r\n";
    text+= "\tor only the selected atoms can be measured, using 'sel'. If 'tot' or 'sel' are\r\n";
    text+= "\tnot specified, the default behavior is 'tot'.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tNsum = <atom count>";

    return text;
}

std::string CCmdCount::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    int num = 0.0;
    int typeOfSelection = 0;

    // Type of selections possible are
    // 0: of total (default), 1: of selected atoms
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text.size() == 0)    typeOfSelection = 0;
    else if(text == "sel")  typeOfSelection = 1;
    else if(text == "tot")  typeOfSelection = 0;
    else
    {
        lastError_ = "Error, type of selection not recognized!";
        return lastError_;
    }

    if(typeOfSelection == 1)
    {
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->isSelected())
                num++;
        }
    }
    else
    {
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            num++;
        }
    }

    fprintf(stdOut_, "\r\n\tNsum = %i\r\n", num);

    return lastError_;
}
