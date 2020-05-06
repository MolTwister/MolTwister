#pragma once
#include "MolTwisterCmd.h"

class CCmdPython : public CCmd
{
public:
    CCmdPython() = delete;
    CCmdPython(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "mtpython"; }
    std::string getTopLevHelpString() const { return std::string("Runs Python script input, with MolTwister spesific functions"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();    
};
