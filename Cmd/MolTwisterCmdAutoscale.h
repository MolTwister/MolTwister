#pragma once
#include "MolTwisterCmd.h"

class CCmdAutoscale : public CCmd
{
public:
    CCmdAutoscale() = delete;
    CCmdAutoscale(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "autoscale"; }
    std::string getTopLevHelpString() const { return std::string("Autoscale the 3D view window"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();    
};
