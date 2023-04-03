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

#include "CmdDefaultAtomProps.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"

std::string CCmdDefaultAtomProps::getCmd()
{
    return "defaultatomprops";
}

std::vector<std::string> CCmdDefaultAtomProps::getCmdLineKeywords()
{
    return { "defaultatomprops", "name", "cpkcolor", "vdwr", "covalentr", "to", "reset" };
}

std::vector<std::string> CCmdDefaultAtomProps::getCmdHelpLines()
{
    return {
                "defaultatomprops cpkcolor name <atom id> to <R> <G> <B>",
                "defaultatomprops vdwr name <atom id> to <r>",
                "defaultatomprops covalentr name <atom id> to <r>",
                "defaultatomprops reset"
           };
}

std::string CCmdDefaultAtomProps::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tChanges the default CPK color, van der Waals radius, and covalent radius, to the value\r\n";
    text+= "\tstated by <R>, <G>, <B>, and <r>. Red, green and blue (i.e., <R>, <G>, <B>) are values\r\n";
    text+= "\tbetween zero and one, while the radius is in units of AA. Note that this only changes\r\n";
    text+= "\tthe values for the currently opened MolTwister instance and will not persist the next\r\n";
    text+= "\ttime MolTwister opens. Hence, if the new default values are used repeatedly, it is\r\n";
    text+= "\trecommended to create a script containing the new values (see 'load' for more information\r\n";
    text+= "\ton scripts).";

    return text;
}

std::string CCmdDefaultAtomProps::execute(std::vector<std::string> arguments)
{
    size_t arg = 0;

    std::string text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "reset")
    {
        state_->defaultAtProp_.resetAtomicPropertiesToDefault();
    }
    else if(text == "cpkcolor")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "name")
        {
            lastError_ = "Syntax Error: Fourth argument should be 'name'!";
            return lastError_;
        }

        std::string atomID = CASCIIUtility::getArg(arguments, arg++);

        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "to")
        {
            lastError_ = "Syntax Error: Sixth argument should be 'to'!";
            return lastError_;
        }

        double R = atof(CASCIIUtility::getArg(arguments, arg++).data());
        double G = atof(CASCIIUtility::getArg(arguments, arg++).data());
        double B = atof(CASCIIUtility::getArg(arguments, arg++).data());
        state_->defaultAtProp_.setCPKColor(atomID, R, G, B);
    }
    else if(text == "vdwr")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "name")
        {
            lastError_ = "Syntax Error: Fourth argument should be 'name'!";
            return lastError_;
        }

        std::string atomID = CASCIIUtility::getArg(arguments, arg++);

        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "to")
        {
            lastError_ = "Syntax Error: Sixth argument should be 'to'!";
            return lastError_;
        }

        double r = atof(CASCIIUtility::getArg(arguments, arg++).data());
        state_->defaultAtProp_.setWDWRadius(atomID, r);
    }
    else if(text == "covalentr")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "name")
        {
            lastError_ = "Syntax Error: Fourth argument should be 'name'!";
            return lastError_;
        }

        std::string atomID = CASCIIUtility::getArg(arguments, arg++);

        text = CASCIIUtility::getArg(arguments, arg++);
        if(text != "to")
        {
            lastError_ = "Syntax Error: Sixth argument should be 'to'!";
            return lastError_;
        }

        double r = atof(CASCIIUtility::getArg(arguments, arg++).data());
        state_->defaultAtProp_.setCovalentRadius(atomID, r);
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'cpkcolor', 'vdwr', or 'covalentr'!";
        return lastError_;
    }

    return lastError_;
}

