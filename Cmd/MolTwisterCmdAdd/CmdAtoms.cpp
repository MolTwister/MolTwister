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

#include "CmdAtoms.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdAtoms::getCmd()
{
    return "atoms";
}

std::vector<std::string> CCmdAtoms::getCmdLineKeywords()
{
    return { "atoms", "sel" };
}

std::vector<std::string> CCmdAtoms::getCmdHelpLines()
{
    return {
                "atoms <ID> sel <dx> <dy> <dz>"
           };
}

std::string CCmdAtoms::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tIn general, ID could for example be C, H, or O (i.e. carbon, hydrogen, or\r\n";
    text+= "\toxygen). It is also possible to give names such as C1 or C2. As long as\r\n";
    text+= "\tthe name contains a C it is recognized as a carbon atom. Similarly, any\r\n";
    text+= "\tname containing O will be recognized as oxygen, etc. \r\n";
    text+= "\r\n";

    text+= "\tWhen atoms are added they attain a new index starting at index 0. The list\r\n";
    text+= "\tof atomic indices are obtained through the 'list' command.\r\n";

    text+= "\r\n";
    text+= "\tThe above command will add one atom of type <ID> a distance <dx> <dy> <dz> away\r\n";
    text+= "\tfrom each selected atom.";

    return text;
}

std::string CCmdAtoms::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CAtom atom;

    text = CASCIIUtility::getArg(arguments, arg++);
    atom.setID(text.data());
    atom.sigma_ = state_->defaultAtProp_.getWDWRadius(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "sel")
    {
        C3DVector   dist;

        text = CASCIIUtility::getArg(arguments, arg++);
        dist.x_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        dist.y_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        dist.z_ = atof(text.data());

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->isSelected() && (state_->currentFrame_ < state_->atoms_[i]->r_.size()))
            {
                atom.r_[0] = state_->atoms_[i]->r_[state_->currentFrame_] + dist;

                state_->addAtom(atom);
            }
        }
    }
    else
    {
        lastError_ = "Syntax Error: Fourth argument should be how to add atoms!";
    }

    return lastError_;
}
