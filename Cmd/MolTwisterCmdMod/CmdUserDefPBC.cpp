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

#include "CmdUserDefPBC.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdUserDefPBC::getCmd()
{
    return "userdefpbc";
}

std::vector<std::string> CCmdUserDefPBC::getCmdLineKeywords()
{
    return { "userdefpbc" };
}

std::vector<std::string> CCmdUserDefPBC::getCmdHelpLines()
{
    return {
                "userdefpbc <x low> <x high> <y low> <y high> <z low> <z high>"
           };
}

std::string CCmdUserDefPBC::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tSets a user defined periodic boundary condition (PBC). Note that the\r\n";
    text+= "\tuser defined PBC will not automatically be applied, but needs to be\r\n";
    text+= "\tactivated through the 'set' command.";

    return text;
}

std::string CCmdUserDefPBC::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    C3DRect pbc;

    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rLow_.x_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rHigh_.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rLow_.y_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rHigh_.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rLow_.z_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    pbc.rHigh_.z_ = atof(text.data());

    if(state_->view3D_) state_->view3D_->setUserPBC(pbc);
    else
    {
        lastError_ = "Could not find 3D View!";
        return lastError_;
    }

    return lastError_;
}

