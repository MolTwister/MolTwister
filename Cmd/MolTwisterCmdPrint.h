#pragma once
#include "MolTwisterCmd.h"

class CCmdPrint : public CCmd
{
public:
    CCmdPrint() = delete;
    CCmdPrint(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "print"; }
    std::string getTopLevHelpString() const { return std::string("Print information (e.g., atom pos. in different file formats)"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseXYZCommand(std::string commandLine, int& arg);
    void parsePDBCommand(std::string commandLine, int& arg);
    void parseMTTCommand(std::string commandLine, int& arg);
    void parseVersionCommand(std::string commandLine, int& arg);
    void parseMixffCommand(std::string commandLine, int& arg);
    bool getMixFFConstituents(std::string commandLine, int& arg, std::vector<std::string>& constituents) const;
};
