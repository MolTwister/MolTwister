#pragma once
#include <vector>
#include "MolTwisterCmd.h"

class CCmdMeasure : public CCmd
{
public:
    CCmdMeasure() = delete;
    CCmdMeasure(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "measure"; }
    std::string getTopLevHelpString() const { return std::string("Perform a measurement"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;
    
protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
