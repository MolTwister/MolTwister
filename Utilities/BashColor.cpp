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

void CBashColor::setColor(EColor foreground, EColor background)
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
    printf("%s", escSeq);
#endif
}

void CBashColor::setSpecial(ESpecial special)
{
#ifndef DEBUG
    char num[3];
    char escSeq[20] = "\033[";

    sprintf(num, "%i", special);
    strcat(escSeq, num);

    strcat(escSeq, "m");
    printf("%s", escSeq);
#endif
}

void CBashColor::clearScreen()
{ 
#ifndef DEBUG
    printf("\033[H\033[J");
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

