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
