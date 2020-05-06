#pragma once
#include "MolTwisterCmd.h"

class CCmdVar : public CCmd
{
public:
    CCmdVar() = delete;
    CCmdVar(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "var"; }
    std::string getTopLevHelpString() const { return std::string("Define a variable"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
};
