#pragma once
#include "MolTwisterCmd.h"

class CCmdMod : public CCmd
{
public:
    CCmdMod() = delete;
    CCmdMod(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "mod"; }
    std::string getTopLevHelpString() const { return std::string("Modify settings, parameters and states"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
