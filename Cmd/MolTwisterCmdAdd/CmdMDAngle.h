#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdMDAngle : public CCmdEntry
{
public:
    CCmdMDAngle() = delete;
    CCmdMDAngle(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdMDAngle() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
