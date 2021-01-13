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

#include "CmdNumCoordinates.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdNumCoordinates::getCmd()
{
    return "numcoordinates";
}

std::vector<std::string> CCmdNumCoordinates::getCmdLineKeywords()
{
    return { "numcoordinates" };
}

std::vector<std::string> CCmdNumCoordinates::getCmdHelpLines()
{
    return {
                "numcoordinates <DCD filename> <record index>"
           };
}

std::string CCmdNumCoordinates::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the number of coordinates in record <record index>\r\n";
    text+= "\tof the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tCoordinate count = <coordinate count>";

    return text;
}

std::string CCmdNumCoordinates::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    int recordIndex;
    std::string text;
    CDCDFile dcdFile;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error opening file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    recordIndex = atoi(text.data());

    dcdFile.gotoRecord(recordIndex);

    fprintf(stdOut_, "\r\n\tCoordinate count = %i\r\n", dcdFile.getNumCoordinatesInRecord());

    return lastError_;
}
