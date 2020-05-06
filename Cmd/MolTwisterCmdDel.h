#pragma once
#include "MolTwisterCmd.h"

class CCmdDel : public CCmd
{
public:
    CCmdDel() = delete;
    CCmdDel(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "del"; }
    std::string getTopLevHelpString() const { return std::string("Delete atoms, bonds, etc."); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();    
};
