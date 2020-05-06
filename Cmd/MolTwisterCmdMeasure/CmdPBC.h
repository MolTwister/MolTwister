#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdPBC : public CCmdEntry
{
public:
    CCmdPBC() = delete;
    CCmdPBC(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdPBC() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
