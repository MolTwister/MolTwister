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

#include "MolTwisterCmdMod.h"

#include "MolTwisterCmdMod/CmdAtomPos.h"
#include "MolTwisterCmdMod/CmdBondLength.h"
#include "MolTwisterCmdMod/CmdAngle.h"
#include "MolTwisterCmdMod/CmdDihedral.h"
#include "MolTwisterCmdMod/CmdCharge.h"
#include "MolTwisterCmdMod/CmdMass.h"
#include "MolTwisterCmdMod/CmdMobillity.h"
#include "MolTwisterCmdMod/CmdSigma.h"
#include "MolTwisterCmdMod/CmdAtomName.h"
#include "MolTwisterCmdMod/CmdResname.h"
#include "MolTwisterCmdMod/CmdUserDefPBC.h"
#include "MolTwisterCmdMod/CmdRotateSel.h"

using namespace MolTwisterCmdMod;

void CCmdMod::onAddKeywords()
{
    addKeyword("mod");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdMod::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdAngle>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdAtomName>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdAtomPos>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdBondLength>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCharge>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDihedral>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMass>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMobility>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdResname>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdRotateSel>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdSigma>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdUserDefPBC>(state_, stdOut_));
}

std::string CCmdMod::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdMod::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdMod::execute(std::string commandLine)
{
    std::string lastError = parser_->executeCmd(commandLine);
    if(!lastError.empty())
    {
        printf("%s\r\n", lastError.data());
        return;
    }

    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1) state_->view3D_->requestUpdate(true);
        else                           state_->view3D_->requestUpdate(false);
    }
}
