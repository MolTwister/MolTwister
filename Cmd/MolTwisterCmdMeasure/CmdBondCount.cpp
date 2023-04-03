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

#include "CmdBondCount.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdBondCount::getCmd()
{
    return "bondcount";
}

std::vector<std::string> CCmdBondCount::getCmdLineKeywords()
{
    return { "bondcount" };
}

std::vector<std::string> CCmdBondCount::getCmdHelpLines()
{
    return {
                "bondcount <atom ID>"
           };
}

std::string CCmdBondCount::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCounts bonds connected to a given atom type, specified by <atom ID> (e.g., H, O, C7), and\r\n";
    text+= "\tpresents a histogram that shows the number of atoms (with atom ID) with 0, 1, 2, .., >=9\r\n";
    text+= "\tbond connections.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. n<=0 n=1 n=2 n=3 n=4 n=5 n=6 n=7 n=8 n>=9\r\n";
    text+= "\t2. <count n<=0> <count n=1> <count n=2> <count n=3> <count n=4> <count n=5> <count n=6> <count n=7> <count n=8> <count n>=9>";

    return text;
}

std::string CCmdBondCount::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    int histogram[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    int bondCount;

    text = CASCIIUtility::getArg(arguments, arg++);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        std::string ID = state_->atoms_[i]->getID();
        if(text != ID) continue;

        bondCount = state_->atoms_[i]->getNumBonds();
        if(bondCount > 9) bondCount = 9;
        if(bondCount < 0) bondCount = 0;

        histogram[bondCount]++;
    }

    printf("\r\n\t%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s%-10s\r\n", "n<=0", "n=1", "n=2", "n=3", "n=4", "n=5", "n=6", "n=7", "n=8", "n>=9");
    printf("\t%-10i%-10i%-10i%-10i%-10i%-10i%-10i%-10i%-10i%-10i\r\n",
           histogram[0], histogram[1], histogram[2], histogram[3], histogram[4], histogram[5], histogram[6], histogram[7], histogram[8], histogram[9]);

    return lastError_;
}
