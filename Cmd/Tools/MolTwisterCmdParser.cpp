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

#include "MolTwisterCmdParser.h"
#include "../../Utilities/ASCIIUtility.h"

std::vector<std::string> CCmdParser::getCmdLineKeywords()
{
    std::vector<std::string> cmdLineKeywords;
    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        std::vector<std::string> cmdEntryKeywords = cmdEntry->getCmdLineKeywords();
        for(auto cmdEntryKeyword : cmdEntryKeywords)
        {
            cmdLineKeywords.emplace_back(cmdEntryKeyword);
        }
    }

    return cmdLineKeywords;
}

std::string CCmdParser::executeCmd(std::string commandLine)
{
    std::string command = CASCIIUtility::getWord(commandLine, 0);
    std::string subCommand = CASCIIUtility::getWord(commandLine, 1);
    if(subCommand.empty()) return std::string("Command ") + command + std::string(" needs an argument!");

    int argIndex = 2;
    std::vector<std::string> argList;
    std::string argText;
    do
    {
        argText = CASCIIUtility::getWord(commandLine, argIndex++);
        if(!argText.empty()) argList.emplace_back(argText);

    } while(!argText.empty());

    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        std::string cmd = cmdEntry->getCmd();
        if(cmd == subCommand)
        {
            std::string errorMsg = cmdEntry->execute(argList);
            if(!errorMsg.empty()) return std::string("Command ") + command + std::string(" ") + subCommand + std::string(" returned error: ") + errorMsg;
            return "";
        }
    }

    return subCommand + std::string(" is not a valid argument for ") + command + std::string("!");
}

std::string CCmdParser::genHelpText(std::string parentCmd)
{
    std::string helpText = std::string("\tOverview of commands of the form '") + parentCmd + " <sub command>':\r\n\r\n";

    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        std::vector<std::string> helpLines = cmdEntry->getCmdHelpLines();
        for(std::string helpLine : helpLines)
        {
            helpText+= std::string("\t") + parentCmd + std::string(" ") + helpLine + std::string("\r\n");
        }
    }

    helpText+= std::string("\r\n\tTo get more information about a <sub command>, type 'help ") + parentCmd + " <sub command>\r\n\r\n";

    return helpText;
}

std::string CCmdParser::genHelpText(std::string parentCmd, std::string subCommand)
{
    std::string helpText = std::string("\tOverview of commands of the form '") + parentCmd + std::string(" ") + subCommand + std::string("':\r\n\r\n");

    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        std::string cmd = cmdEntry->getCmd();
        if(cmd == subCommand)
        {
            std::vector<std::string> helpLines = cmdEntry->getCmdHelpLines();
            for(std::string helpLine : helpLines)
            {
                helpText+= std::string("\t") + parentCmd + std::string(" ") + helpLine + std::string("\r\n");
            }

            helpText+= std::string("\r\n") + cmdEntry->getCmdFreetextHelp() + std::string("\r\n\r\n");

            return helpText;
        }
    }

    return std::string("\t") + subCommand + std::string(" is not a valid sub command argument for ") + parentCmd + std::string("!\r\n\r\n");
}

std::shared_ptr<std::vector<std::string>> CCmdParser::getListOfCmdEntryCommands()
{
    auto cmdStringList = std::make_shared<std::vector<std::string>>();

    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        cmdStringList->emplace_back(cmdEntry->getCmd());
    }

    return cmdStringList;
}

void CCmdParser::redirectOutput(FILE* stdOut)
{
    for(std::shared_ptr<CCmdEntry> cmdEntry : cmdEntryList_)
    {
        cmdEntry->redirectOutput(stdOut);
    }
}
