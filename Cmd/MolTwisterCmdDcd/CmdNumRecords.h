#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdNumRecords : public CCmdEntry
{
public:
    CCmdNumRecords() = delete;
    CCmdNumRecords(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdNumRecords() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
