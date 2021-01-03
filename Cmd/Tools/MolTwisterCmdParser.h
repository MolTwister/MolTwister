#pragma once
#include "MolTwisterCmdEntry.h"

class CCmdParser
{
public:
    CCmdParser() = default;

public:
    void registerCmd(std::shared_ptr<CCmdEntry> cmdEntry) { cmdEntryList_.emplace_back(cmdEntry); }
    std::vector<std::string> getCmdLineKeywords();
    std::string executeCmd(std::string commandLine);
    std::string genHelpText(std::string parentCmd);
    std::string genHelpText(std::string parentCmd, std::string subCommand);
    std::shared_ptr<std::vector<std::string>> getListOfCmdEntryCommands();
    void purge() { cmdEntryList_.clear(); }
    void redirectOutput(FILE* stdOut = stdout);

private:
    std::vector<std::shared_ptr<CCmdEntry>> cmdEntryList_;
};
