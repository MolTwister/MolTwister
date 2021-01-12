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

#include "CmdNumRecords.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdNumRecords::getCmd()
{
    return "numrecords";
}

std::vector<std::string> CCmdNumRecords::getCmdLineKeywords()
{
    return { "numrecords" };
}

std::vector<std::string> CCmdNumRecords::getCmdHelpLines()
{
    return {
                "numrecords <DCD filename>"
           };
}

std::string CCmdNumRecords::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the number of recoords within the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tNum records = <record count>";

    return text;
}

std::string CCmdNumRecords::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CDCDFile dcdFile;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error opening file ") + text + std::string("!");
        return lastError_;
    }

    fprintf(stdOut_, "\r\n\tNum records = %i\r\n", dcdFile.getNumRecords());

    return lastError_;
}
