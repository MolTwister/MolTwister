#pragma once
#include "MolTwisterCmd.h"

class CCmdAdd : public CCmd
{
public:
    CCmdAdd() = delete;
    CCmdAdd(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "add"; }
    std::string getTopLevHelpString() const { return std::string("Add atom to system"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
