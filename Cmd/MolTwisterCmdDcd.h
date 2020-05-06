#pragma once
#include "MolTwisterCmd.h"

class CCmdDcd : public CCmd
{
public:
    CCmdDcd() = delete;
    CCmdDcd(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "dcd"; }
    std::string getTopLevHelpString() const { return std::string("Perform a DCD file operation"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();
};
