#pragma once
#include "MolTwisterCmd.h"

class CCmdVarlist : public CCmd
{
public:
    CCmdVarlist() = delete;
    CCmdVarlist(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "varlist"; }
    std::string getTopLevHelpString() const { return std::string("List all variables defined by 'var'"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();    
};
