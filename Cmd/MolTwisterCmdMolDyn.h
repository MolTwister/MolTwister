#pragma once
#include "MolTwisterCmd.h"
#include "MolTwisterCmdMolDyn/Config/MolDynConfig.h"

class CCmdMolDyn : public CCmd
{
public:
    CCmdMolDyn() = delete;
    CCmdMolDyn(CMolTwisterState* state) : CCmd(state) { state_ = state; }

public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "moldyn"; }
    std::string getTopLevHelpString() const { return std::string("Invoke a molecular dynamics simulation"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();

private:
    static CMolDynConfig molDynConfig_;
};
