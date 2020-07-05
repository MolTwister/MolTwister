#pragma once
#include "../Tools/MolTwisterCmdEntry.h"
#include "Config/MolDynConfig.h"

class CCmdCfg : public CCmdEntry
{
public:
    CCmdCfg() = delete;
    CCmdCfg(CMolTwisterState* state, FILE* stdOut, CMolDynConfig* molDynConfig) : CCmdEntry(state, stdOut) { molDynConfig_ = molDynConfig; }
    virtual ~CCmdCfg() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void parseGetCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseSetCommand(const std::vector<std::string>& arguments, size_t& arg);
    void parseResettodefaultsCommand(const std::vector<std::string>& arguments, size_t& arg);

private:
    std::string lastError_;
    CMolDynConfig* molDynConfig_;
};
