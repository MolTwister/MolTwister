#pragma once
#include "MolTwisterCmd.h"

class CCmdHoomdblue : public CCmd
{
public:
    CCmdHoomdblue() = delete;
    CCmdHoomdblue(CMolTwisterState* state) : CCmd(state) { state_ = state; }

public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "hoomdblue"; }
    std::string getTopLevHelpString() const { return std::string("Operations related to the HOOMD-blue SW package"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseGenffCommand(std::string commandLine, int& arg);
    void parseGendataCommand(std::string commandLine, int& arg);
    void parseGenrunCommand(std::string commandLine, int& arg);
    void parseAtmtoljCommand(std::string commandLine, int& arg);
    void parseFstoljCommand(std::string commandLine, int& arg);
    void parseKtoljCommand(std::string commandLine, int& arg);
};
