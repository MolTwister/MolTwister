#include <iostream>
#include <errno.h>
#include <dirent.h>
#include <unistd.h>
#include "MolTwisterCmdCd.h"

void CCmdCd::onAddKeywords()
{
    addKeyword("cd");
}

std::string CCmdCd::getHelpString() const
{ 
    std::string text;
    
    text+= "\tUsage: cd <destination directory>\r\n";
    text+= "\r\n";
    text+= "\tChanges directory to <destination directory>. It is possible to specify\r\n";
    text+= "\tdirectory names containing spaces using the '\\' character. For example,\r\n";
    text+= "\t'cd /Users/John\\ Doe/' would change to the directory '/Users/John Doe'\r\n";
    text+= "\tTo change directory to MolTwister shortcut N, type 'cd [N]' and hit enter.\r\n";
    text+= "\tShortcuts can be edited directly by accessing 'MolTwister.shortcuts', \r\n";
    text+= "\tlocated under the home directoy of your computer. The contents of this\r\n";
    text+= "\tfile is structured as follows:\r\n";
    text+= "\r\n";
    text+= "\t\t#Default\r\n";
    text+= "\t\t<default directory at startup>\r\n";
    text+= "\t\t#n\r\n";
    text+= "\t\t<shortcut 1>\r\n";
    text+= "\t\t     .\r\n";
    text+= "\t\t     .\r\n";
    text+= "\t\t     .\r\n";
    text+= "\t\t<shortcut n>\r\n";
    text+= "\r\n";
    text+= "\twhere n is the number of shortcuts available in the file.";
    
    return text;
}

void CCmdCd::execute(std::string commandLine)
{
    std::string text;
    
    // Process MolTwister version of 'cd'
    text = getDirWord(commandLine, 1);
    
    if(text.find('[') != std::string::npos)
    {
        CASCIIUtility::removeWhiteSpace(text, "[]");
        int shortcutIndex = (int)atoi(text.data()) - 1;
        
        if((shortcutIndex < 0) || (shortcutIndex >= state_->shortcutDirs_.size()))
        {
            printf("Error: Invalid shortcut index!");
        }
        else
        {
            if(state_)
            {
                if(chdir(state_->shortcutDirs_[shortcutIndex].data()) == -1)
                {
                    printf("cd %s failed - %s", state_->shortcutDirs_[shortcutIndex].data(), strerror(errno));
                }
            }
        }
    }
    
    else
    {
        if(chdir(text.data()) == -1)
        {
            printf("cd failed - %s", strerror(errno));
        }
    }
}

std::string CCmdCd::getDirWord(std::string line, int wordIndex, const char* whiteSpaceChars) const
{
    size_t lineLen = line.size();
    bool inAWord = false;
    char prevChar;
    int currIndex;
    
    if(CASCIIUtility::isWhiteSpace(line[0], whiteSpaceChars))
    {
        prevChar = ' ';
        inAWord = false;
        currIndex = -1;
    }
    else
    {
        prevChar = line[0];
        inAWord = true;
        currIndex = 0;
    }
    
    std::string word;
    for(int i=0; i<lineLen; i++)
    {
        if(inAWord)
        {
            if(CASCIIUtility::isWhiteSpace(line[i], whiteSpaceChars) && (prevChar != '\\'))
            {
                inAWord = false;
            }
            else
            {
                if(wordIndex == currIndex)
                {
                    prevChar = line[i];
                    
                    if(prevChar != '\\') word+= line[i];
                }
            }
        }
        else
        {
            if(!CASCIIUtility::isWhiteSpace(line[i], whiteSpaceChars))
            {
                inAWord = true;
                currIndex++;
                
                if(wordIndex == currIndex)
                {
                    word+= line[i];
                }
            }            
        }
    }

    return word;
}
