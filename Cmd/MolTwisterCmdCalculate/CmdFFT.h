#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdFFT : public CCmdEntry
{
public:
    CCmdFFT() = delete;
    CCmdFFT(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdFFT() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
