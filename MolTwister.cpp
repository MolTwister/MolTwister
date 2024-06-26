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
#include <pthread.h>
#include <math.h>
#include <stdlib.h>
#include <filesystem>
#include <readline/readline.h>
#include <readline/history.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <string>
#include <vector>
#include "Utilities/BashColor.h"
#include "MolTwisterCommandPool.h"
#include "MolTwister.h"
#include "3DView/MolTwister3DView.h"

std::vector<std::shared_ptr<CCmd>> CMolTwister::cmdList_;


CMolTwister::CMolTwister()
    : tutorialPool_({ "./moltwister-tutorials/", "/usr/local/bin/moltwister-tutorials/" })
{
    startDir_ = dirWorking;
}

CMolTwister::~CMolTwister()
{
}

void CMolTwister::run(C3DView* view3D, const std::string& hostIP, const std::string& hostPort)
{
    pthread_t threadHandle;

    state_.view3D_ = view3D;
    hostIP_ = hostIP;
    hostPort_ = hostPort;
    
    if(pthread_create(&threadHandle, nullptr, threadRun, (void*)this))
    {
        printf("Error: Could not create main thread!\r\n");
    }
}

std::string CMolTwister::genTerminalStartupInfo() const
{
    std::string introText = CBashColor::clearScreen(false);

    introText+= CBashColor::setSpecial(CBashColor::specBright, false);
    introText+= "----------------------------------------------------------------\r\n";
    introText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
    introText+= CBashColor::setColor(CBashColor::colYellow, CBashColor::colNone, false);
    introText+= CASCIIUtility::sprintf("                        MolTwister V%s                       \r\n", MOLTWISTER_VER);
    introText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
    introText+= CBashColor::setSpecial(CBashColor::specBright, false);
    introText+= "----------------------------------------------------------------\r\n";
    introText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
    introText+= CBashColor::setColor(CBashColor::colYellow, CBashColor::colNone, false);
    introText+= " Copyright (C) 2023 Richard Olsen.\r\n";
    introText+= "\r\n";
    introText+= " This program is free software: you can redistribute it and/or\r\n";
    introText+= " modify it under the terms of the GNU General Public License as\r\n";
    introText+= " published by the Free Software Foundation, either version 3 of\r\n";
    introText+= " the License, or (at your option) any later version.\r\n";
    introText+= "\r\n";
    introText+= " This program is distributed in the hope that it will be useful,\r\n";
    introText+= " but WITHOUT ANY WARRANTY; without even the implied warranty of\r\n";
    introText+= " MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\r\n";
    introText+= " GNU General Public License for more details.\r\n";
    introText+= "\r\n";
    introText+= " You should have received a copy of the GNU General Public\r\n";
    introText+= " License along with this program. If not, see\r\n";
    introText+= "               <https://www.gnu.org/licenses/>.\r\n";
    introText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
    introText+= "\r\n";
    introText+= " Write 'help' to list available commands, write 'help <command>'\r\n";
    introText+= " to obtain help for a specific command, and type 'exit' to exit\r\n";
    introText+= " the program. \r\n";
    introText+= "\r\n";
    introText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);

    return introText;
}

std::string CMolTwister::genHelp() const
{
    std::string helpText = "\r\n";

    for(int i=0; i<(int)cmdList_.size(); i++)
    {
        helpText+= CBashColor::setColor(CBashColor::colGreen, CBashColor::colNone, false);
        helpText+= CBashColor::setSpecial(CBashColor::specBright, false);
        helpText+= CASCIIUtility::sprintf("\t%-15s", cmdList_[i]->getCmd().data());
        helpText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
        helpText+= CASCIIUtility::sprintf(" - %s\r\n", cmdList_[i]->getTopLevHelpString().data());
    }

    helpText+= CASCIIUtility::sprintf("\r\n");
    for(int i=0; i<(int)tutorialPool_.getCount(); i++)
    {
        const std::string tutorialString = std::string("tutorial") + std::to_string(i+1);
        const std::string tutorialDescription = tutorialPool_.getDocumentHeader(i);

        helpText+= CBashColor::setColor(CBashColor::colMagenta, CBashColor::colNone, false);
        helpText+= CBashColor::setSpecial(CBashColor::specBright, false);
        helpText+= CASCIIUtility::sprintf("\t%-15s", tutorialString.data());
        helpText+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
        helpText+= CASCIIUtility::sprintf(" - %s\r\n", tutorialDescription.data());
    }

    helpText+= CASCIIUtility::sprintf("\r\n\tAny commands not recognized by MolTwister will be\r\n");
    helpText+= CASCIIUtility::sprintf("\tredirected to the shell (e.g. ssh, vim, pico). Write\r\n");
    helpText+= CASCIIUtility::sprintf("\t'help __MarkDown__' to create a Markdown help document.\r\n");

    helpText+= CASCIIUtility::sprintf("\r\n\tIn the 3D View you can use the mouse to:\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Rotate the system by holding down the left mouse\r\n");
    helpText+= CASCIIUtility::sprintf("\t     button inside the view while moving the mouse\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Zoom in and out by holding down the middle mouse\r\n");
    helpText+= CASCIIUtility::sprintf("\t     button (or mouse wheel) inside the view while moving\r\n");
    helpText+= CASCIIUtility::sprintf("\t     the mouse\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Pan the system by holding down the right mouse button\r\n");
    helpText+= CASCIIUtility::sprintf("\t     inside the view while moving the mouse\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Select or deselect atoms by holding down shift, while\r\n");
    helpText+= CASCIIUtility::sprintf("\t     clicking the left mouse button\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Deselect all atoms by holding down shift, while clicking\r\n");
    helpText+= CASCIIUtility::sprintf("\t     the right mouse button\r\n");
    helpText+= CASCIIUtility::sprintf("\r\n");
    helpText+= CASCIIUtility::sprintf("\tHotkeys when 3D view is the active window:\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Alt + f: Fullscreen (use Esc to exit fullscreen mode)\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Alt + o: Switch between orthographic and perspective mode\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Alt + a: Switch on or off axis\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Alt + i: Switch on or off Coulomb energy iso-surfaces.\r\n");
    helpText+= CASCIIUtility::sprintf("\t              Should only be used on smaller systems, with a\r\n");
    helpText+= CASCIIUtility::sprintf("\t              few molecules, due to the slow nature of the\r\n");
    helpText+= CASCIIUtility::sprintf("\t              applied algorithms.\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Alt + p: Switch on or off visibility of bonds stretching\r\n");
    helpText+= CASCIIUtility::sprintf("\t              across periodic boundaries (PBC)\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * ->:      Go forward one frame or 10, 100, 1000, with Shift,\r\n");
    helpText+= CASCIIUtility::sprintf("\t              Alt, Shift+Alt, respectively\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * <-:      Go back one frame or 10, 100, 1000, with Shift,\r\n");
    helpText+= CASCIIUtility::sprintf("\t              Alt, Shift+Alt, respectively\r\n");
    helpText+= CASCIIUtility::sprintf("\r\n");
    helpText+= CASCIIUtility::sprintf("\tNotes on units and conventions:\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * The Angstrom unit is denoted as AA\r\n");
    helpText+= CASCIIUtility::sprintf("\r\n");
    helpText+= CASCIIUtility::sprintf("\tTutorials:\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * To view a tutorial, execute the appropirate 'tutorial#' command\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * The tutorial is written in Markdown, which is human readable text\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * You can store the Markdown text to file by 'tutorial# > filename.md'\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * Use a Markdown compiler to convert to HTML or PDF.\r\n");
    helpText+= CASCIIUtility::sprintf("\t     E.g., Visual Studio Code and a Markdown plugin\r\n");
    helpText+= CASCIIUtility::sprintf("\t   * The tutorials are also included in the Markdown file created by\r\n");
    helpText+= CASCIIUtility::sprintf("\t     the 'help __MarkDown__' command\r\n");

    return helpText;
}

std::string CMolTwister::genHelpForCommand(const std::string& cmdString, const std::string& commandLine) const
{
    std::string helpOnCommand;

    int cmdIndex = -1;

    for(int i=0; i<(int)cmdList_.size(); i++)
    {
        if(cmdList_[i]->checkCmd(cmdString.data()))
        {
            cmdIndex = i;
            break;
        }
    }

    if(cmdIndex >= 0)
    {
        std::string subCmdString = CASCIIUtility::getWord(commandLine, 2);
        CASCIIUtility::removeWhiteSpace(subCmdString, " \t\r\n");

        helpOnCommand+= CBashColor::setColor(CBashColor::colGreen, CBashColor::colNone, false);
        helpOnCommand+= CBashColor::setSpecial(CBashColor::specBright, false);
        helpOnCommand+= CASCIIUtility::sprintf("\r\n\t%s", cmdList_[cmdIndex]->getCmd().data());
        helpOnCommand+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
        if(!subCmdString.empty())
        {
            helpOnCommand+= CBashColor::setColor(CBashColor::colCyan, CBashColor::colNone, false);
            helpOnCommand+= CASCIIUtility::sprintf(" %s", subCmdString.data());
            helpOnCommand+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
        }
        helpOnCommand+= CASCIIUtility::sprintf(" - %s\r\n", cmdList_[cmdIndex]->getTopLevHelpString().data());

        helpOnCommand+= CBashColor::setSpecial(CBashColor::specBright, false);
        helpOnCommand+= CASCIIUtility::sprintf("\t---------------------------------------------------------------------------\r\n\r\n");
        helpOnCommand+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);

        if(subCmdString.empty())
            helpOnCommand+= CASCIIUtility::sprintf("%s\r\n", cmdList_[cmdIndex]->getHelpString().data());
        else
            helpOnCommand+= CASCIIUtility::sprintf("%s\r\n", cmdList_[cmdIndex]->getHelpString(subCmdString).data());


        helpOnCommand+= CBashColor::setSpecial(CBashColor::specBright, false);
        helpOnCommand+= CASCIIUtility::sprintf("\t---------------------------------------------------------------------------\r\n");
        helpOnCommand+= CBashColor::setColor(CBashColor::colNone, CBashColor::colNone, false);
    }
    else
    {
        helpOnCommand+= CASCIIUtility::sprintf("\r\n\tCould not find help for %s!\r\n", cmdString.data());
    }

    return helpOnCommand;
}

std::string CMolTwister::genMarkdownHelpDoc() const
{
    std::string markdownDoc;

    // List commands
    markdownDoc+= CASCIIUtility::sprintf("# List of commands\r\n\r\n");
    for(int i=0; i<(int)cmdList_.size(); i++)
    {
        std::string markdownText;

        markdownDoc+= CASCIIUtility::sprintf("## %s\r\n\r\n", cmdList_[i]->getCmd().data());
        markdownText = cmdList_[i]->getTopLevHelpString();
        markdownText = CASCIIUtility::trimFromLeft(markdownText, "\t", 1);
        markdownDoc+= CASCIIUtility::sprintf("*%s*\r\n\r\n", markdownText.data());

        markdownText = cmdList_[i]->getHelpString();
        markdownText = CASCIIUtility::createMarkDownCodeBlock(markdownText, 4, true);
        markdownDoc+= CASCIIUtility::sprintf("%s\r\n", markdownText.data());

        std::shared_ptr<std::vector<std::string>> subCmdList = cmdList_[i]->getListOfSubCommands();
        if(subCmdList)
        {
            for(int j=0; j<(int)subCmdList->size(); j++)
            {
                markdownDoc+= CASCIIUtility::sprintf("### %s\r\n\r\n", (*subCmdList)[j].data());

                markdownText = cmdList_[i]->getHelpString((*subCmdList)[j]);
                markdownText = CASCIIUtility::createMarkDownCodeBlock(markdownText, 4, true);
                markdownDoc+= CASCIIUtility::sprintf("%s\r\n", markdownText.data());
            }
        }
    }

    // List tutorials
    for(int i=0; i<tutorialPool_.getCount(); i++)
    {
        markdownDoc+= CASCIIUtility::sprintf("\r\n%s\r\n", tutorialPool_.getDocument(i).data());
    }

    return markdownDoc;
}

void* CMolTwister::threadRun(void* arg)
{
    CMolTwister* mtPtr = (CMolTwister*)arg;
    
    if(!mtPtr->_run())
    {
        printf("MolTwister failed!!!\r\n\r\n");
    }
    
    return nullptr;
}

void CMolTwister::initReadline() const
{
    // Allow conditional parsing of the ~/.inputrc file.
    rl_readline_name = (char*)("MolTwister");
    
    // Tell the completer that we want a crack first.
    rl_attempted_completion_function = commandCompletionFunc;
}

char** CMolTwister::commandCompletionFunc(const char* text, int, int)
{
    char** matchesPtr;
    
    matchesPtr = rl_completion_matches(text, commandCompletionGenerator);
    
    return (matchesPtr);
}

char* CMolTwister::duplicateString(const char* str)
{
    char* stringRet;
    
    stringRet = (char*)malloc(strlen(str) + 1);
    strcpy(stringRet, str);

    return stringRet;
}

char* CMolTwister::commandCompletionGenerator(const char* text, int state)
{
    static int keyIndex, len;
    std::string name;
    
    // If this is a new word to complete, initialize now. This includes
    // saving the length of TEXT for efficiency, and initializing the index
    // variable to 0.
    if(!state)
    {
        keyIndex = 0;
        len = (int)strlen(text);
    }
    
    // Return the next name which partially matches from the command list.
    for(int j=keyIndex; j<CCmd::getNumKeywords(); j++)
    {
        name = CCmd::getKeyword(j);
        keyIndex++;

        if(strncmp(name.data(), text, len) == 0)
            return duplicateString(name.data());
    }

    // If no names matched, then return nullptr.
    return nullptr;
}

std::string CMolTwister::readLine() const
{    
    // Retrieve leaf dir. to show in prompt
    std::string absDir = CFileUtility::getCWD();
    std::string leafDir = CFileUtility::extractLeafDir(absDir);
    
    // Show prompt and read user input
    std::string commandLine;

    std::string prompt = "MolTwister:";
    prompt+= leafDir;
    prompt+= "> ";
    commandLine = readline(prompt.data());
    if(commandLine.size() > 0) add_history(commandLine.data());

    return commandLine;
}

bool CMolTwister::_run()
{
    std::function<std::pair<FILE*, std::string>(const int&, const std::string&, std::function<void(FILE*)>)> handlePiping =
            [](const int& pipeSymbIndex, const std::string& commandLine, std::function<void(FILE* fileStdOut)> onRedirectOutput)
    {
        // Check if '> <file>' (i.e. pipe) was requested. If so
        // redirect all fprintf() output to <file> and keep all printf()
        // statements directed to stdout
        std::string fileName;
        FILE* fileStdOut = stdout;
        if(pipeSymbIndex != -1)
        {
            std::string commmandLineCpy = commandLine;

            fileName = commmandLineCpy.substr(pipeSymbIndex+1, std::string::npos);
            CASCIIUtility::removeWhiteSpace(fileName);
            fileStdOut = fopen(fileName.data(), "w+");

            if(fileStdOut)  onRedirectOutput(fileStdOut);
            else            printf("Error: could not create file!");
        }

        return std::pair<FILE*, std::string>(fileStdOut, fileName);
    };

    std::function<void(FILE*, const std::string&, const bool&)> endPiping = [](FILE* fileStdOut, const std::string& fileName, const bool& skipTabRemove)
    {
        // Close file, if redirection was requested, then
        // remove unwanted '\t' (i.e. tab) from that file
        if(fileStdOut != stdout)
        {
            fclose(fileStdOut);
            if(!skipTabRemove) CFileUtility::removeTabsFromFile(fileName);
        }
    };

    std::string command;
    std::string commandLine;
    bool quit = false;

    printf("%s", genTerminalStartupInfo().c_str());

    cmdList_.clear();
    CMolTwisterCommandPool::generateCmdList(&state_, cmdList_);
    
    for(int i=0; i<(int)cmdList_.size(); i++) cmdList_[i]->init();
    
    readShortcutDirs();
    initReadline();

    do
    {        
        // Read from command line
        CBashColor::setSpecial(CBashColor::specBright);
        commandLine = readLine();
        CBashColor::setColor();
        if(commandLine.size() == 0) continue;

        // Retrieve command
        commandLine = CASCIIUtility::removeCRLF(commandLine);
        command = CASCIIUtility::getWord(commandLine, 0);
        
        // Parse command
        if(command == "exit")
        {
            if(state_.view3D_) state_.view3D_->requestQuit();
            quit = true;
            printf("\r\nExiting MolTwister...\r\n");
        }
        
        else if(command == "help")
        {            
            std::string cmdString = CASCIIUtility::getWord(commandLine, 1);
            CASCIIUtility::removeWhiteSpace(cmdString, " \t\r\n");
            
            if(cmdString.size() == 0)
            {
                printf("%s", genHelp().c_str());
            }
            else if(cmdString == "__MarkDown__")
            {
                std::string mkdFileName = "list_of_commands.md";
                FILE* mkdFile = fopen(mkdFileName.data(), "w");

                if(mkdFile)
                {
                    fprintf(mkdFile, "%s", genMarkdownHelpDoc().c_str());

                    fclose(mkdFile);
                    printf("\r\n\tCreated file %s!\r\n", mkdFileName.data());
                }
                else
                {
                    printf("\r\n\tError! could not create file %s\r\n", mkdFileName.data());
                }
            }
            else
            {
                printf("%s", genHelpForCommand(cmdString, commandLine).c_str());
            }
        }
        
        else
        {
            bool foundCmd = false;
            int pipeSymbIndex = CASCIIUtility::findString("> ", commandLine);

            for(int i=0; i<(int)cmdList_.size(); i++)
            {
                if(cmdList_[i] && cmdList_[i]->checkCmd(command.data()))
                {
                    foundCmd = true;

                    // Do piping of command to desired location, if requested
                    auto [fileStdOut, fileName] = handlePiping(pipeSymbIndex, commandLine, [this, i](FILE* fileStdOut) { cmdList_[i]->redirectOutput(fileStdOut); });

                    // Execute command
                    cmdList_[i]->execute(commandLine);
                    cmdList_[i]->redirectOutput(stdout);
                    bool skipTabRemove = cmdList_[i]->checkSkipTabRemoveReqFlagAndReset();

                    // End piping, if it was started
                    endPiping(fileStdOut, fileName, skipTabRemove);
                }                
            }

            for(int i=0; i<tutorialPool_.getCount(); i++)
            {
                std::string tutorialCmd = std::string("tutorial") + std::to_string(i+1);

                if(command == tutorialCmd)
                {
                    foundCmd = true;

                    // Do piping of command to desired location, if requested
                    auto [fileStdOut, fileName] = handlePiping(pipeSymbIndex, commandLine, [this, i](FILE* fileStdOut){ fprintf(fileStdOut, "%s", tutorialPool_.getDocument(i).data()); });

                    // Print to stdout if this is the destination
                    if(fileStdOut == stdout) printf("\r\n%s\r\n", CASCIIUtility::addTabsToDocument(tutorialPool_.getDocument(i)).data());

                    // End piping, if it was started
                    endPiping(fileStdOut, fileName, true);
                    break;
                }
            }

            if(!foundCmd)
            {
                if(system(nullptr))
                {
                    printf("\r\n");
                    if(system(commandLine.data()) != -1)
                    {
                        foundCmd = true;
                    }
                }
                
                if(!foundCmd)
                {
                    printf("Unknown command %s, type 'help' to list available commands!", command.data());
                }
            }
        }
        
        printf("\r\n");

    } while(!quit);

    cmdList_.clear();
    return true;
}

void CMolTwister::readShortcutDirs()
{
    struct passwd* pw;
    std::string line = "~";
    std::string currWorkingDir;
    FILE* fileHandle;
    bool isLastLineInFile;

    currWorkingDir = CFileUtility::getCWD();
    
    pw = getpwuid(getuid());
    if(pw == nullptr) return;
    const char* szHomedir = pw->pw_dir;
    chdir(szHomedir);

    fileHandle = fopen("MolTwister.shortcuts", "r");
    if(fileHandle)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, isLastLineInFile);
            if(line.find("#") != std::string::npos)
            {
                CASCIIUtility::removeWhiteSpace(line, " \t#\r\n");
                if(line.find("Default") != std::string::npos)
                {
                    line = CFileUtility::readLine(fileHandle, isLastLineInFile);
                    if(startDir_ == dirInitFile)
                    {
                        if(chdir(line.data()) == -1)
                        {
                            printf("Error while reading MolTwister.shortcuts: Could not locate directory %s!", line.data());
                        }
                    }
                }
                else
                {
                    int iNumDirectories = atoi(line.data());
                    
                    for(int i=0; i<iNumDirectories; i++)
                    {
                        line = CFileUtility::readLine(fileHandle, isLastLineInFile);
                        CASCIIUtility::removeWhiteSpace(line, "\r\n");
                        state_.shortcutDirs_.emplace_back(line);
                    }

                    isLastLineInFile = true;
                }
            }
            else
            {
                printf("Error while reading MolTwister.shortcuts: First line must be either #Default or #<num>!");
                break;
            }

        } while(!isLastLineInFile);
        
        fclose(fileHandle);
    }

    if(startDir_ == dirWorking) chdir(currWorkingDir.data());
}
