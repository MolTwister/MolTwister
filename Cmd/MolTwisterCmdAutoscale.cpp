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

#include <iostream>
#include "MolTwisterCmdAutoscale.h"

void CCmdAutoscale::onAddKeywords()
{
    addKeyword("autoscale");
}

std::string CCmdAutoscale::getHelpString() const
{ 
    std::string text;
    
    text+= "\tUsage: autoscale\r\n";
    text+= "\r\n";
    text+= "\tAutoscale the 3D view by searching for the boundary of all atoms inside the\r\n";
    text+= "\t3D view and subsequently adjusting the viewing frustum to these boundaries.\r\n";
    text+= "\tThe camera is placed along the x-axis and will point towards the center of\r\n";
    text+= "\tthe located atomic boundaries.\r\n";
    
    return text;
}

void CCmdAutoscale::execute(std::string)
{
    if(!state_) return;

    if(state_->view3D_) state_->view3D_->requestUpdate(true);
}
