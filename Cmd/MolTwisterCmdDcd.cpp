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

#include "MolTwisterCmdDcd.h"

#include "MolTwisterCmdDcd/CmdReadRecord.h"
#include "MolTwisterCmdDcd/CmdNumRecords.h"
#include "MolTwisterCmdDcd/CmdReadCoordinate.h"
#include "MolTwisterCmdDcd/CmdNumCoordinates.h"
#include "MolTwisterCmdDcd/CmdHeader.h"
#include "MolTwisterCmdDcd/CmdUnwrap.h"
#include "MolTwisterCmdDcd/CmdAtomicUnwrap.h"
#include "MolTwisterCmdDcd/CmdWrap.h"
#include "MolTwisterCmdDcd/CmdFromXTC.h"

void CCmdDcd::onAddKeywords()
{
    addKeyword("dcd");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdDcd::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdAtomicUnwrap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdHeader>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdNumCoordinates>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdNumRecords>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdReadCoordinate>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdReadRecord>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdUnwrap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdWrap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdFromXTC>(state_, stdOut_));
}

std::string CCmdDcd::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdDcd::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdDcd::execute(std::string commandLine)
{
    std::string lastError = parser_->executeCmd(commandLine);
    if(!lastError.empty())
    {
        printf("%s\r\n", lastError.data());
        return;
    }
}
