#pragma once
#include "MolTwisterCmd.h"

class CCmdLammps : public CCmd
{    
public:
    CCmdLammps() = delete;
    CCmdLammps(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "lammps"; }
    std::string getTopLevHelpString() const { return std::string("Operations related to the LAMMPS SW package"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseGenffCommand(std::string commandLine, int& arg);
    void parseGendataCommand(std::string commandLine, int& arg);
};
