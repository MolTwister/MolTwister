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

#include <iostream>
#include <vector>
#include "MolTwisterCmdSet.h"

void CCmdSet::onAddKeywords()
{
    addKeyword("set");
    addKeyword("projection");
    addKeyword("perspective");
    addKeyword("ortho");
    addKeyword("bondacrosspbc");
    addKeyword("fullscreen");
    addKeyword("on");
    addKeyword("off");
    addKeyword("userdefpbc");
    addKeyword("redrawlimit");
    addKeyword("fog");
}

std::string CCmdSet::getHelpString() const
{
    std::string text;

    text+= "\tUsage: set projection <projection>\r\n";
    text+= "\t       set fullscreen <fullscreen>\r\n";
    text+= "\t       set userdefpbc <userdefpbc>\r\n";
    text+= "\t       set bondacrosspbc <bondacrosspbc>\r\n";
    text+= "\t       set redrawlimit <num atoms>\r\n";
    text+= "\t       set fog <fog>\r\n";
    text+= "\r\n";
    text+= "\tThis command will set various properties of MolTwister or its loaded systems.\r\n";
    text+= "\r\n";
    text+= "\tThe <projection> can be either 'perspective' or 'ortho'.\r\n";
    text+= "\r\n";
    text+= "\tThe <fullscreen> parameter can be either 'on' or 'off'.\r\n";
    text+= "\r\n";
    text+= "\tThe <userdefpbc> parameter can be either 'on' or 'off'. If 'on', the user defined\r\n";
    text+= "\tperiodic boundary conditions (PBC) set by the 'mod userdefpbc' command is applied.\r\n";
    text+= "\tIf 'off', the PBC will be collected elsewhere.\r\n";
    text+= "\r\n";
    text+= "\tThe <bondacrosspbc> parameter can be either 'on' or 'off'. If 'on', bonds acrsss\r\n";
    text+= "\tPBCs will be visible in the 3D view.\r\n";
    text+= "\r\n";
    text+= "\tThe 'redrawlimit' sub-command will specify that, during rotation / panning / scaling\r\n";
    text+= "\tof the scene, scenes with more than <num atoms> atoms will only be displayed with the\r\n";
    text+= "\taxis system. All atoms will be redrawn once the corresponding mouse button is released.\r\n";
    text+= "\r\n";
    text+= "\tThe <fog> can be either 'on' or 'off'.";

    return text;
}

void CCmdSet::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "projection")
    {
        parseProjectionCommand(commandLine, arg);
    }

    else if(text == "fullscreen")
    {
        parseFullscreenCommand(commandLine, arg);
    }

    else if(text == "userdefpbc")
    {
        parseUserdefpbcCommand(commandLine, arg);
    }

    else if(text == "bondacrosspbc")
    {
        parseBondacrosspbcCommand(commandLine, arg);
    }

    else if(text == "redrawlimit")
    {
        parseRedrawlimitCommand(commandLine, arg);
    }

    else if(text == "fog")
    {
        parseFogCommand(commandLine, arg);
    }

    else
    {
        printf("Syntax Error: Second argument should specify what to set!");
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1)  state_->view3D_->requestUpdate(true);
        else                            state_->view3D_->requestUpdate(false);
    }
}

void CCmdSet::parseProjectionCommand(std::string commandLine, int& arg)
{
    std::string text;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "perspective")
    {
        if(state_->view3D_)
        {
            state_->view3D_->setOrthoView(false);
        }
    }
    else if(text == "ortho")
    {
        if(state_->view3D_)
        {
            state_->view3D_->setOrthoView(true);
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'perspective' or 'ortho'!");
    }
}

void CCmdSet::parseBondacrosspbcCommand(std::string commandLine, int& arg)
{
    std::string text;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "on")
    {
        if(state_->view3D_)
        {
            state_->view3D_->enableViewOfBondsAcrossPBC(true);
        }
    }
    else if(text == "off")
    {
        if(state_->view3D_)
        {
            state_->view3D_->enableViewOfBondsAcrossPBC(false);
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'on' or 'off'!");
    }
}

void CCmdSet::parseFullscreenCommand(std::string commandLine, int& arg)
{
    std::string text;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "on")
    {
        if(state_->view3D_)
        {
            state_->view3D_->requestFullScreen(true);
        }
    }
    else if(text == "off")
    {
        if(state_->view3D_)
        {
            state_->view3D_->requestFullScreen(false);
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'on' or 'off'!");
    }
}

void CCmdSet::parseUserdefpbcCommand(std::string commandLine, int& arg)
{
    std::string text;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "on")
    {
        if(state_->view3D_)
        {
            state_->view3D_->enableUserPBC(true);
        }
    }
    else if(text == "off")
    {
        if(state_->view3D_)
        {
            state_->view3D_->enableUserPBC(false);
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'on' or 'off'!");
    }
}

void CCmdSet::parseRedrawlimitCommand(std::string commandLine, int& arg)
{
    C3DView::setRedrawLimit(atoi(CASCIIUtility::getWord(commandLine, arg++).data()));
}

void CCmdSet::parseFogCommand(std::string commandLine, int& arg)
{
    std::string text;

    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "on")
    {
        if(state_->view3D_)
        {
            state_->view3D_->setFog(true);
        }
    }
    else if(text == "off")
    {
        if(state_->view3D_)
        {
            state_->view3D_->setFog(false);
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'on' or 'off'!");
    }
}
