#include <iostream>
#include <sys/ioctl.h>
#include <stdio.h>
#include <errno.h>
#include <dirent.h>
#include <string>
#include <unistd.h>
#include "Utilities/BashColor.h"
#include "MolTwisterCmdLs.h"

void CCmdLs::onAddKeywords()
{
    addKeyword("lsm");
}

std::string CCmdLs::getHelpString() const
{ 
    std::string  text;
    
    text+= "\tUsage: lsm\r\n";
    text+= "\r\n";
    text+= "\tShows the contents of the current directory in the same way as the 'ls'\r\n";
    text+= "\tcommand in Linux, but will also show the directory shortcuts that are\r\n";
    text+= "\tavailable to the MolTwister 'cd' command. Go to the help pages for the\r\n";
    text+= "\tMolTwister 'cd' command for a further explanation on shortcuts. However, \r\n";
    text+= "\tthe 'lsm' command does not accept any of theswitches available for the\r\n";
    text+= "\tLinux 'ls' command.";
    
    return text;
}

void CCmdLs::execute(std::string)
{
    struct winsize terminalSize;
    std::string text;
    int numColumns = 3;
    int lineIndex = 0;
    
    ioctl(0, TIOCGWINSZ, &terminalSize);
    if(terminalSize.ws_col > 0)
        numColumns = (terminalSize.ws_col - strlen("\t")) / 50;
    
    text = CFileUtility::getCWD();
    if(text.empty())
        printf("Could not retrieve current working directory!");
    else
    {
        DIR* dirHandle = opendir(text.data());
        dirent* dirEntry;
        
        if(dirHandle)
        {
            fprintf(stdOut_, "\r\n");
            
            do
            {
                dirEntry = readdir(dirHandle);
                if(dirEntry)
                {
                    if(lineIndex % numColumns == 0)
                        fprintf(stdOut_, "\t");
                    
                    if(dirEntry->d_type == DT_DIR)
                    {
                        CBashColor::setSpecial(CBashColor::specBright);
                        CBashColor::setColor(CBashColor::colGreen);
                        fprintf(stdOut_, "[DIR] %-44s", dirEntry->d_name);
                        CBashColor::setColor();
                    }
                    else
                    {
                        fprintf(stdOut_, "%-50s", dirEntry->d_name);
                    }

                    if(lineIndex % numColumns == (numColumns-1))
                        fprintf(stdOut_, "\r\n");
                    
                    lineIndex++;
                }
                
            } while(dirEntry);

            
            for(int i=0; i<state_->shortcutDirs_.size(); i++)
            {
                if(i == 0) fprintf(stdOut_, "\r\n\r\n");
                CBashColor::setColor(CBashColor::colBlue);
                fprintf(stdOut_, "\t[%i] - %s\r\n", i+1, state_->shortcutDirs_[i].data());
                CBashColor::setColor();
            }

            fprintf(stdOut_, "\r\n");
        }
        else
            printf("Could not open handle to current directory!");
        
        closedir(dirHandle);
    }
}
