#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdUserDefPBC : public CCmdEntry
{
public:
    CCmdUserDefPBC() = delete;
    CCmdUserDefPBC(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdUserDefPBC() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
