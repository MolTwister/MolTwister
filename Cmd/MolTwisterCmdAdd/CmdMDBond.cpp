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

#include "CmdMDBond.h"
#include "../../Utilities/ASCIIUtility.h"

#define DEFAULT_MAX_BOND 1.75

std::string CCmdMDBond::getCmd()
{
    return "mdbond";
}

std::vector<std::string> CCmdMDBond::getCmdLineKeywords()
{
    return { "mdbond" };
}

std::vector<std::string> CCmdMDBond::getCmdHelpLines()
{
    return {
                "mdbond <ID1> <ID2> <FF-type> [all_r_less_than <radius>, mol_r_less_than <radius>, only_visible_bonds, only_14_bonds] <parameters for given FF-type>"
           };
}

std::string CCmdMDBond::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tID1 and ID2 identifiy the atom types where a bond interaction is to be\r\n";
    text+= "\tdefined (e.g., H, O, C, O5). The force field type, FF-type, is idefined by\r\n";
    text+= "\ta string. The possible strings are listed below, together with how 'parameters\r\n";
    text+= "\tfor given FF-type' is defined for each string.\r\n";
    text+= "\r\n";
    text+= "\tPossible bond definitions to apply:\r\n";
    text+= "\t* all_r_less_than <r> - yields bonds even between atoms not close enough to define a molecule (if the r-criteria is satisfied)\r\n";
    text+= "\t* mol_r_less_than <r> - yields bonds only between atoms close enough to define a molecule (if the r-criteria is satisfied)\r\n";
    text+= "\t* only_visible_bonds - yields bonds only where bonds are visible in the 3D view\r\n";
    text+= "\t* only_14_bonds - yields bonds only between 1-4 interactions\r\n";

    text+= "\r\n";

    text+= "\tPossible <FF-type> <parameters for given FF-type> to apply:\r\n";

    int regTypes = state_->mdFFBondList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        CMDFFBond* ffType = state_->mdFFBondList_.getRegisteredFFType(i);
        if(ffType)
        {
            std::string ffTypeString = ffType->getFFType();
            std::vector<std::string> paramsHelpLines = ffType->getCmdHelpLines();

            for(std::string paramsHelp : paramsHelpLines)
            {
                text+= std::string("\t* ") + ffTypeString + std::string(" ") + paramsHelp + std::string("\r\n");
            }
        }
    }

    return text;
}

std::string CCmdMDBond::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CMDFFBond::EDetCrit detCrit;
    std::shared_ptr<CMDFFBond> bondEntry;
    std::string stringAtom1, stringAtom2, type, text;
    double detCritR0;
    int regTypes;

    stringAtom1 = CASCIIUtility::getArg(arguments, arg++);
    stringAtom2 = CASCIIUtility::getArg(arguments, arg++);
    type = CASCIIUtility::getArg(arguments, arg++);

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "all_r_less_than")
    {
        detCrit = CMDFFBond::critAllRLessThan_R0;
        text = CASCIIUtility::getArg(arguments, arg++);
        detCritR0 = atof(text.data());
    }
    else if(text == "mol_r_less_than")
    {
        detCrit = CMDFFBond::critMolRLessThan_R0;
        text = CASCIIUtility::getArg(arguments, arg++);
        detCritR0 = atof(text.data());
    }
    else if(text == "only_visible_bonds")
    {
        detCrit = CMDFFBond::critOnlyVisibleBonds;
        detCritR0 = 0.0;
    }
    else if(text == "only_14_bonds")
    {
        detCrit = CMDFFBond::critOnly14Bonds;
        detCritR0 = 0.0;
    }
    else
    {
        detCrit = CMDFFBond::critAllRLessThan_R0;
        detCritR0 = DEFAULT_MAX_BOND;
        arg--;
    }

    regTypes = state_->mdFFBondList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        if(type == state_->mdFFBondList_.getRegisteredFFType(i)->getFFType())
        {
            bondEntry = state_->mdFFBondList_.getRegisteredFFType(i)->createCopy();
            break;
        }
    }

    if(bondEntry)
    {
        bondEntry->setAtomsToBond(stringAtom1, stringAtom2);
        bondEntry->setBondDetectionCriteria(detCrit, detCritR0);

        arguments.erase(arguments.begin(), arguments.begin() + arg);
        bondEntry->parse(arguments);
        state_->mdFFBondList_.add(*bondEntry);
    }
    else
    {
        lastError_ = std::string("Error: Did not recognize bond type ") + type + std::string("!");
        return lastError_;
    }

    return lastError_;
}
