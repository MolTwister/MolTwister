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

#include "CmdMDNonBonded.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdMDNonBonded::getCmd()
{
    return "mdnonbonded";
}

std::vector<std::string> CCmdMDNonBonded::getCmdLineKeywords()
{
    return { "mdnonbonded" };
}

std::vector<std::string> CCmdMDNonBonded::getCmdHelpLines()
{
    return {
                "mdnonbonded <ID1> <ID2> <FF-type> <parameters for given FF-type>"
           };
}

std::string CCmdMDNonBonded::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tID1 and ID2 identifiy the atom types where a non-bonded interaction is to\r\n";
    text+= "\tbe defined (e.g., H, O, C, O5). The force field type, FF-type, is idefined by\r\n";
    text+= "\ta string. The possible strings are listed below, together with how 'parameters\r\n";
    text+= "\tfor given FF-type' is defined for each string.\r\n";
    text+= "\r\n";

    text+= "\tPossible <FF-type> <parameters for given FF-type> to apply:\r\n";

    int regTypes = state_->mdFFNonBondedList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        CMDFFNonBonded* ffType = state_->mdFFNonBondedList_.getRegisteredFFType(i);
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

std::string CCmdMDNonBonded::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::shared_ptr<CMDFFNonBonded> nonBondedEntry;
    std::string stringAtom1, stringAtom2, text;

    stringAtom1 = CASCIIUtility::getArg(arguments, arg++);
    stringAtom2 = CASCIIUtility::getArg(arguments, arg++);
    text = CASCIIUtility::getArg(arguments, arg++);

    int regTypes = state_->mdFFNonBondedList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        if(text == state_->mdFFNonBondedList_.getRegisteredFFType(i)->getFFType())
        {
            nonBondedEntry = state_->mdFFNonBondedList_.getRegisteredFFType(i)->createCopy();
            break;
        }
    }

    if(nonBondedEntry)
    {
        nonBondedEntry->setAtomsToBond(stringAtom1, stringAtom2);

        std::vector<std::string> ffArguments = arguments;
        ffArguments.erase(ffArguments.begin(), ffArguments.begin() + arg);
        if(!nonBondedEntry->parse(ffArguments))
        {
            lastError_ = std::string("Error parsing: ") + CASCIIUtility::argsToString(arguments) + std::string("!");
            return lastError_;
        }
        state_->mdFFNonBondedList_.add(*nonBondedEntry);
    }
    else
    {
        lastError_ = std::string("Error: Did not recognize non-bonded type ") + text + std::string("!");
        return lastError_;
    }

    return lastError_;
}
