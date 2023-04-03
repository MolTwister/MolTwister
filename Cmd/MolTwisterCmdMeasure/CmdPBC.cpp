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

#include "CmdPBC.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdPBC::getCmd()
{
    return "pbc";
}

std::vector<std::string> CCmdPBC::getCmdLineKeywords()
{
    return { "pbc" };
}

std::vector<std::string> CCmdPBC::getCmdHelpLines()
{
    return {
                "pbc [<frame index>]"
           };
}

std::string CCmdPBC::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures the current periodic bounary conditions of frame <frame index>.\r\n";
    text+= "\tThe frame index is, by default, the current frame index (i.e., the visible\r\n";
    text+= "\tframe), in case <frame index> is omitted.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. x = [<x low>, <x high>]\r\n";
    text+= "\t2. x = [<y low>, <y high>]\r\n";
    text+= "\t3. x = [<z low>, <z high>]\r\n";
    text+= "\tIf the PBC is user defined, this will be notified by the message 'User\r\n";
    text+= "\tdefined PBC!', in addition to the above.";

    return text;
}

std::string CCmdPBC::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    C3DRect pbc;
    std::string text;

    text = CASCIIUtility::getArg(arguments, arg++);

    if(text.size() == 0)
        pbc = state_->view3D_->calcPBC(state_->currentFrame_);
    else
    {
        pbc = state_->view3D_->calcPBC(atoi(text.data()));
    }

    fprintf(stdOut_, "\r\n\tx = [%.4f, %.4f]\r\n", pbc.rLow_.x_, pbc.rHigh_.x_);
    fprintf(stdOut_, "\ty = [%.4f, %.4f]\r\n", pbc.rLow_.y_, pbc.rHigh_.y_);
    fprintf(stdOut_, "\tz = [%.4f, %.4f]\r\n", pbc.rLow_.z_, pbc.rHigh_.z_);

    if(state_->view3D_->isUserPBCEnabled())
        fprintf(stdOut_, "\r\n\tUser defined PBC!\r\n");

    return lastError_;
}
