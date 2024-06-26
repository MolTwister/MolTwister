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

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "BashColor.h"

#define CLI_HIDE_CUR
#define CLI_SHOW_CUR
#define CLI_SAVE_CUR_POS
#define CLI_REST_CUR_POS
#define CLI_CLR_LAST_LINE
#define CLI_MV_CUR_UP
#define CLI_MV_CUR_DN

std::string CBashColor::setColor(EColor foreground, EColor background, bool printToStdOut)
{
#ifndef DEBUG
    char num[3];
    char escSeq[20] = "\033[";
    
    if (!foreground && !background) strcat(escSeq, "0");     // Reset colors if no params
    
    if(foreground)
    {
        sprintf(num, "%i", 29 + foreground);
        strcat(escSeq, num);
        
        if(background) strcat(escSeq, ";");
    }
    
    if(background)
    {
        sprintf(num, "%i", 39 + background);
        strcat(escSeq, num);
    }
    
    strcat(escSeq, "m");
    if(printToStdOut) printf("%s", escSeq);
    return escSeq;
#else
    return "";
#endif
}

std::string CBashColor::setSpecial(ESpecial special, bool printToStdOut)
{
#ifndef DEBUG
    char num[3];
    char escSeq[20] = "\033[";

    sprintf(num, "%i", special);
    strcat(escSeq, num);

    strcat(escSeq, "m");
    if(printToStdOut) printf("%s", escSeq);
    return escSeq;
#else
    return "";
#endif
}

std::string CBashColor::clearScreen(bool printToStdOut)
{ 
#ifndef DEBUG
    std::string escSeq = "\033[H\033[J";
    if(printToStdOut) printf("%s", escSeq.c_str());
    return escSeq;
#else
    return "";
#endif
}



void CBashControl::saveCursorPos()
{
#ifndef DEBUG
    printf("\033[s");
#endif
}

void CBashControl::restoreCursorPos()
{
#ifndef DEBUG
    printf("\033[u");
#endif
}

void CBashControl::hideCursor()
{
#ifndef DEBUG
    printf("\033[?25l");
#endif
}

void CBashControl::showCursor()
{
#ifndef DEBUG
    printf("\033[?25h");
#endif
}

void CBashControl::clearLastLine()
{
#ifndef DEBUG
    printf("\033[A\033[2K");
#endif
}

void CBashControl::moveCursorUp()
{
#ifndef DEBUG
    printf("\033[A");
#endif
}

void CBashControl::moveCursorDown()
{
#ifndef DEBUG
    printf("\033[B");
#endif
}

void CBashControl::moveCursorTo(int x, int y)
{
#ifndef DEBUG
    printf("\033[%i;%iH", y, x);
#endif
}
