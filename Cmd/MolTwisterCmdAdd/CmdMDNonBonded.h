#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdMDNonBonded : public CCmdEntry
{
public:
    CCmdMDNonBonded() = delete;
    CCmdMDNonBonded(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdMDNonBonded() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
