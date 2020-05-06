#pragma once
#include <vector>
#include "MolTwisterCmd.h"

class CCmdCalculate : public CCmd
{
public:
    CCmdCalculate() = delete;
    CCmdCalculate(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    virtual void execute(std::string commandLine);
    virtual std::string getCmd() const { return "calculate"; }
    std::string getTopLevHelpString() const { return std::string("Perform a calculation"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
