#pragma once
#include "MolTwisterCmd.h"

class CCmdCopy : public CCmd
{
public:
    CCmdCopy() = delete;
    CCmdCopy(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "copy"; }
    std::string getTopLevHelpString() const { return std::string("Copy atoms, molecules, etc."); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
};
