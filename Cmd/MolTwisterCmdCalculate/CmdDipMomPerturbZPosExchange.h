#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdDipMomPerturbZPosExchange : public CCmdEntry
{
public:
    CCmdDipMomPerturbZPosExchange() = delete;
    CCmdDipMomPerturbZPosExchange(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdDipMomPerturbZPosExchange() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
