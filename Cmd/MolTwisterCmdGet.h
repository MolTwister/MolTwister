#pragma once
#include "MolTwisterCmd.h"

class CCmdGet : public CCmd
{
public:
    CCmdGet() = delete;
    CCmdGet(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "get"; }
    std::string getTopLevHelpString() const { return std::string("Get various properties of MolTwister and the loaded systems"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
    
private:
    void parseAtomtypesCommand(std::string commandLine, int& arg);
    void parseMdinconsistencyCommand(std::string commandLine, int& arg);
    void parseBondinfoCommand(std::string commandLine, int& arg);
    void parseUserdefpbcCommand(std::string commandLine, int& arg);
    void parseGpuinfoCommand(std::string commandLine, int& arg);
};
