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

#include <iostream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <errno.h>
#include <dirent.h>
#include <string>
#include <unistd.h>
#include "Utilities/BashColor.h"
#include "MolTwisterCmdLl.h"

void CCmdLl::onAddKeywords()
{
    addKeyword("ll");
}

std::string CCmdLl::getHelpString() const
{ 
    std::string  text;
    
    text+= "\tUsage: ll\r\n";
    text+= "\r\n";
    text+= "\tThis is a shortcut for the 'ls -lha' command (only applicable on Unix based systems).";
    
    return text;
}

void CCmdLl::execute(std::string)
{
    printf("\r\n");
    system("ls -lha");
}
