#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdCenter : public CCmdEntry
{
public:
    CCmdCenter() = delete;
    CCmdCenter(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdCenter() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
