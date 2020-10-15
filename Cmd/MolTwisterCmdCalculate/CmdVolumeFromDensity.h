#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdVolumeFromDensity : public CCmdEntry
{
public:
    CCmdVolumeFromDensity() = delete;
    CCmdVolumeFromDensity(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdVolumeFromDensity() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
