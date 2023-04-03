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

#include "CmdRotateSel.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"

std::string CCmdRotateSel::getCmd()
{
    return "rotatesel";
}

std::vector<std::string> CCmdRotateSel::getCmdLineKeywords()
{
    return { "rotatesel" };
}

std::vector<std::string> CCmdRotateSel::getCmdHelpLines()
{
    return {
                "rotatesel <angle> <pos. of rotation vector> <rotation vector>"
           };
}

std::string CCmdRotateSel::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tRotates the current selection <angle> degrees around the vector <rotation vector>,\r\n";
    text+= "\twhich is positioned at <pos. of rotation vector>. The position and vector are formatted\r\n";
    text+= "\tas three numbers, <x> <y> <z>, separated by space.";

    return text;
}

std::string CCmdRotateSel::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;

    // Get angle to rotate
    text = CASCIIUtility::getArg(arguments, arg++);
    double angle = atof(text.data());
    angle*= M_PI / 180.0;

    // Get pos of vector to rotate around
    C3DVector pos;
    text = CASCIIUtility::getArg(arguments, arg++);
    pos.x_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    pos.y_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    pos.z_ = atof(text.data());

    // Get vector to rotate around
    C3DVector vec;
    text = CASCIIUtility::getArg(arguments, arg++);
    vec.x_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vec.y_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vec.z_ = atof(text.data());

    // Rotate selection
    C3DBasis basis;
    CMolecularTools molTools(state_, stdOut_);
    basis.generateCartessianBasisAt(pos, vec);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]->isSelected())
        {
            C3DVector r = state_->atoms_[i]->r_[state_->currentFrame_];
            C3DVector r_new = molTools.rotatePosAroundBasisW(pos, basis, r, angle);
            state_->atoms_[i]->r_[state_->currentFrame_] = r_new;
        }
    }

    return lastError_;
}

