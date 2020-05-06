#pragma once
#include <string>
#include "MolTwisterCmd.h"

class CCmdCd : public CCmd
{
public:
    CCmdCd() = delete;
    CCmdCd(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "cd"; }
    std::string getTopLevHelpString() const { return std::string("Change directory (MolTwister implementation)"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();
    
private:
    std::string getDirWord(std::string line, int wordIndex, const char* whiteSpaceChars=nullptr) const;
};
