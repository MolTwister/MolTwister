#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdQBal : public CCmdEntry
{
public:
    CCmdQBal() = delete;
    CCmdQBal(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdQBal() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void findAtomTypesInSel(std::vector<std::string>& types) const;
    void findInstancesOfAtomType(std::string type, bool includeSelOnly, std::vector<CAtom*>& atoms) const;
    double sgn(double x) const { if(x < 0.0) return -1.0; if(x > 0.0) return 1.0; return 0.0; }

private:
    std::string lastError_;
};
