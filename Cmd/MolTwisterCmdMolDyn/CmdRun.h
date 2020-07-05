#pragma once
#include "../Tools/MolTwisterCmdEntry.h"
#include "Config/MolDynConfig.h"

class CCmdRun : public CCmdEntry
{
public:
    CCmdRun() = delete;
    CCmdRun(CMolTwisterState* state, FILE* stdOut, CMolDynConfig* molDynConfig) : CCmdEntry(state, stdOut) { molDynConfig_ = molDynConfig; }
    virtual ~CCmdRun() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
    CMolDynConfig* molDynConfig_;
};
