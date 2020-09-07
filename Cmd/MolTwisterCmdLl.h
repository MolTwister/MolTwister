#pragma once
#include "MolTwisterCmd.h"

class CCmdLl : public CCmd
{
public:
    CCmdLl() = delete;
    CCmdLl(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "ll"; }
    std::string getTopLevHelpString() const { return std::string("Shortcut for the Unix 'ls -lha' command"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
};
