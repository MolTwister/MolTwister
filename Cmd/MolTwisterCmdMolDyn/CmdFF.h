#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdFF : public CCmdEntry
{
public:
    CCmdFF() = delete;
    CCmdFF(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdFF() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void parseBondforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseAngleforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseDihedralforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseNonbondforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg);

private:
    std::string lastError_;
};
