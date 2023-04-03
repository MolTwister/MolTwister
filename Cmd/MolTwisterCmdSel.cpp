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
#include <vector>
#include "Utilities/ConditionalOnXYZ.h"
#include "MolTwisterCmdSel.h"

void CCmdSel::onAddKeywords()
{
    addKeyword("sel");
    addKeyword("atom");
    addKeyword("atomname");
    addKeyword("within");
    addKeyword("all");
    addKeyword("none");
}

std::string CCmdSel::getHelpString() const
{
    std::string text;

    text+= "\tUsage: sel atom <atom index>\r\n";
    text+= "\t       sel atomname <atom ID (e.g., O, H, C7)> [within <condition>]\r\n";
    text+= "\t       sel all [within <domain>]\r\n";
    text+= "\t       sel none\r\n";
    text+= "\r\n";
    text+= "\tThis command will create an atomic selection, where the selected atoms will be highlighted\r\n";
    text+= "\tin the 3D view. The first time a selection is done, the selection is done. The second time\r\n";
    text+= "\tthe same selection is done, the selection is unselected, except for 'sel all' and 'sel none'.\r\n";
    text+= "\r\n";
    text+= CConditionalOnXYZ::getHelpText("<domain>", false);

    return text;
}

void CCmdSel::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "atom")
    {
        parseAtomCommand(commandLine, arg);
    }

    else if(text == "atomname")
    {
        parseAtomnameCommand(commandLine, arg);
    }
    
    else if(text == "all")
    {
        parseAllCommand(commandLine, arg);
    }

    else if(text == "none")
    {
        parseNoneCommand(commandLine, arg);
    }
    
    else
    {
        printf("Syntax Error: Second argument should specify what to select!");
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1)  state_->view3D_->requestUpdate(true);
        else                                state_->view3D_->requestUpdate(false);
    }
}

void CCmdSel::parseAtomCommand(std::string commandLine, int& arg)
{
    std::string text;
    int atomIndex;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    atomIndex = atoi(text.data());
    
    if(atomIndex < state_->atoms_.size())
    {
        if(state_->atoms_[atomIndex]->isSelected())
            state_->atoms_[atomIndex]->select(false);
        else
            state_->atoms_[atomIndex]->select(true);
    }
    else
    {
        printf("Error: could not find requested atom!");
    }
}

void CCmdSel::parseAtomnameCommand(std::string commandLine, int& arg)
{
    CConditionalOnXYZ cond;
    std::string name, text;
    bool isWithin = false;
    
    name = CASCIIUtility::getWord(commandLine, arg++);
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "within")
    {
        // Domain specifier of the form x<5&x>8|y<=4&y>3...
        // without any space in-between
        text = CASCIIUtility::getWord(commandLine, arg++);
        cond.parse(text);
        isWithin = true;
    }
    
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        std::string ID = state_->atoms_[i]->getID();
        
        if(ID == name)
        {
            if(!isWithin)
            {
                if(state_->atoms_[i]->isSelected())
                    state_->atoms_[i]->select(false);
                else
                    state_->atoms_[i]->select(true);
            }
            else
            {
                if(state_->currentFrame_ >= state_->atoms_[i]->r_.size()) continue;
                C3DVector v = state_->atoms_[i]->r_[state_->currentFrame_];
                if(cond.isFulfilled(v.x_, v.y_, v.z_))
                {
                    if(state_->atoms_[i]->isSelected())
                        state_->atoms_[i]->select(false);
                    else
                        state_->atoms_[i]->select(true);
                }
            }
        }
    }
}

void CCmdSel::parseAllCommand(std::string commandLine, int& arg)
{
    CConditionalOnXYZ cond;
    std::string text;
    bool within = false;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "within")
    {
        // Domain specifier of the form x<5&x>8|y<=4&y>3...
        // without any space in-between
        text = CASCIIUtility::getWord(commandLine, arg++);
        cond.parse(text);
        within = true;
    }

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(!within)
        {
            state_->atoms_[i]->select(true);
        }
        else
        {
            if(state_->currentFrame_ >= state_->atoms_[i]->r_.size()) continue;
            C3DVector v = state_->atoms_[i]->r_[state_->currentFrame_];
            if(cond.isFulfilled(v.x_, v.y_, v.z_))
            {
                state_->atoms_[i]->select(true);
            }
        }
    }
}

void CCmdSel::parseNoneCommand(std::string, int&)
{
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        state_->atoms_[i]->select(false);
    }
}
