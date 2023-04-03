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

#include "MolTwisterCommandPool.h"

#include "Cmd/MolTwisterCmdCd.h"
#include "Cmd/MolTwisterCmdLs.h"
#include "Cmd/MolTwisterCmdLl.h"
#include "Cmd/MolTwisterCmdGenBonds.h"
#include "Cmd/MolTwisterCmdList.h"
#include "Cmd/MolTwisterCmdAutoscale.h"
#include "Cmd/MolTwisterCmdDel.h"
#include "Cmd/MolTwisterCmdCopy.h"
#include "Cmd/MolTwisterCmdAdd.h"
#include "Cmd/MolTwisterCmdVar.h"
#include "Cmd/MolTwisterCmdVarlist.h"
#include "Cmd/MolTwisterCmdMod.h"
#include "Cmd/MolTwisterCmdSet.h"
#include "Cmd/MolTwisterCmdGet.h"
#include "Cmd/MolTwisterCmdSel.h"
#include "Cmd/MolTwisterCmdPrint.h"
#include "Cmd/MolTwisterCmdMeasure.h"
#include "Cmd/MolTwisterCmdCalculate.h"
#include "Cmd/MolTwisterCmdGauss9.h"
#include "Cmd/MolTwisterCmdLammps.h"
#include "Cmd/MolTwisterCmdHoomdBlue.h"
#include "Cmd/MolTwisterCmdMmol.h"
#include "Cmd/MolTwisterCmdDcd.h"
#include "Cmd/MolTwisterCmdPython.h"
#include "Cmd/MolTwisterCmdLoad.h"
#include "Cmd/MolTwisterCmdMolDyn.h"

void CMolTwisterCommandPool::generateCmdList(CMolTwisterState* mtState, std::vector<std::shared_ptr<CCmd>>& cmdList)
{
    cmdList.emplace_back(std::make_shared<CCmdCd>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdLs>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdLl>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdLoad>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdGenBonds>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdList>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdAutoscale>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdDel>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdCopy>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdAdd>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdVar>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdVarlist>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdMod>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdSet>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdGet>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdSel>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdPrint>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdMeasure>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdCalculate>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdGauss9>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdLammps>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdHoomdblue>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdMmol>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdDcd>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdPython>(mtState));
    cmdList.emplace_back(std::make_shared<CCmdMolDyn>(mtState));

    for(size_t i=0; i<cmdList.size(); i++) cmdList[i]->init();
}
