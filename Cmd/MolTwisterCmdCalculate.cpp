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

#include <stdio.h>
#include "MolTwisterCmdCalculate.h"

#include "MolTwisterCmdCalculate/CmdQBal.h"
#include "MolTwisterCmdCalculate/CmdDensityProfile.h"
#include "MolTwisterCmdCalculate/CmdVACF.h"
#include "MolTwisterCmdCalculate/CmdFFT.h"
#include "MolTwisterCmdCalculate/CmdVDOS.h"
#include "MolTwisterCmdCalculate/CmdMSD.h"
#include "MolTwisterCmdCalculate/CmdPairCorrelation.h"
#include "MolTwisterCmdCalculate/CmdCOM.h"
#include "MolTwisterCmdCalculate/CmdPotenEnergyMap.h"
#include "MolTwisterCmdCalculate/CmdForceBetween.h"
#include "MolTwisterCmdCalculate/CmdEnergyBetween.h"
#include "MolTwisterCmdCalculate/CmdDensityMap.h"
#include "MolTwisterCmdCalculate/CmdLoading.h"
#include "MolTwisterCmdCalculate/CmdEnergyOfTranslation.h"
#include "MolTwisterCmdCalculate/CmdHBondCount.h"
#include "MolTwisterCmdCalculate/CmdDihedralDistr.h"
#include "MolTwisterCmdCalculate/CmdDihedralDistrCOM.h"
#include "MolTwisterCmdCalculate/CmdDistProbabilityCOM.h"
#include "MolTwisterCmdCalculate/CmdDipMomProfile.h"
#include "MolTwisterCmdCalculate/CmdDipMomPerturbCharge.h"
#include "MolTwisterCmdCalculate/CmdDipMomPerturbZPosExchange.h"
#include "MolTwisterCmdCalculate/CmdDipMomPerturbZPos.h"
#include "MolTwisterCmdCalculate/CmdDipMom.h"
#include "MolTwisterCmdCalculate/CmdVolumeFromDensity.h"

void CCmdCalculate::onAddKeywords()
{
    addKeyword("calculate");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdCalculate::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdCOM>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDensityMap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDensityProfile>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDihedralDistr>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDihedralDistrCOM>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDipMom>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDipMomPerturbCharge>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDipMomPerturbZPos>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDipMomPerturbZPosExchange>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDipMomProfile>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDistProbabilityCOM>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdEnergyBetween>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdEnergyOfTranslation>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdFFT>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdForceBetween>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdHBondCount>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdLoading>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMSD>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdPairCorrelation>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdPotenEnergyMap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdQBal>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdVACF>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdVDOS>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdVolumeFromDensity>(state_, stdOut_));
}

std::string CCmdCalculate::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdCalculate::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdCalculate::execute(std::string commandLine)
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
