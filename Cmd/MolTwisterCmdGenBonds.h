#pragma once
#include "MolTwisterCmd.h"

class CCmdGenBonds : public CCmd
{
public:
    CCmdGenBonds() = delete;
    CCmdGenBonds(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "genbonds"; }
    std::string getTopLevHelpString() const { return std::string("Generate bonds between atoms"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();    
};
