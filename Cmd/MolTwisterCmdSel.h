#pragma once
#include "MolTwisterCmd.h"

class CCmdSel : public CCmd
{
public:
    CCmdSel() = delete;
    CCmdSel(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "sel"; }
    std::string getTopLevHelpString() const { return std::string("Create a selection of atoms"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseAtomCommand(std::string commandLine, int& arg);
    void parseAtomnameCommand(std::string commandLine, int& arg);
    void parseAllCommand(std::string commandLine, int& arg);
    void parseNoneCommand(std::string commandLine, int& arg);
};
