#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdMDBond : public CCmdEntry
{
public:
    CCmdMDBond() = delete;
    CCmdMDBond(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdMDBond() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
