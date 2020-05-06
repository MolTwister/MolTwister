#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdNumCoordinates : public CCmdEntry
{
public:
    CCmdNumCoordinates() = delete;
    CCmdNumCoordinates(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdNumCoordinates() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    std::string lastError_;
};
