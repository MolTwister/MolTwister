#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdAtoms : public CCmdEntry
{
public:
    CCmdAtoms() = delete;
    CCmdAtoms(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdAtoms() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
