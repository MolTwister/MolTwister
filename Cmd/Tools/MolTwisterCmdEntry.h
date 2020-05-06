#pragma once
#include "../../MolTwisterState.h"
#include <string>
#include <vector>

class CCmdEntry
{
public:
    CCmdEntry() = delete;
    CCmdEntry(CMolTwisterState* state, FILE* stdOut) { state_ = state; stdOut_ = stdOut; }
    virtual ~CCmdEntry() = default;

public:
    virtual std::string getCmd() = 0;
    virtual std::vector<std::string> getCmdLineKeywords() = 0;
    virtual std::vector<std::string> getCmdHelpLines() = 0;
    virtual std::string getCmdFreetextHelp() = 0;
    virtual std::string execute(std::vector<std::string> arguments) = 0;

protected:
    CMolTwisterState* state_;
    FILE* stdOut_;
};
