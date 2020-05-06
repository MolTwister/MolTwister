#pragma once
#include "MolTwisterCmd.h"

class CCmdSet : public CCmd
{
public:
    CCmdSet() = delete;
    CCmdSet(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "set"; }
    std::string getTopLevHelpString() const { return std::string("Set various properties of MolTwister and the loaded systems"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();

private:
    void parseProjectionCommand(std::string commandLine, int& arg);
    void parseFullscreenCommand(std::string commandLine, int& arg);
    void parseUserdefpbcCommand(std::string commandLine, int& arg);
    void parseBondacrosspbcCommand(std::string commandLine, int& arg);
    void parseRedrawlimitCommand(std::string commandLine, int& arg);
    void parseFogCommand(std::string commandLine, int& arg);
};
