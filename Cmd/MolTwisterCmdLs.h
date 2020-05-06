#pragma once
#include "MolTwisterCmd.h"

class CCmdLs : public CCmd
{
public:
    CCmdLs() = delete;
    CCmdLs(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "lsm"; }
    std::string getTopLevHelpString() const { return std::string("Show current directory contents"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
};
