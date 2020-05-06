#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

class CCmdAtom : public CCmdEntry
{
public:
    CCmdAtom() = delete;
    CCmdAtom(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
    virtual ~CCmdAtom() = default;

public:
    std::string getCmd();
    std::vector<std::string> getCmdLineKeywords();
    std::vector<std::string> getCmdHelpLines();
    std::string getCmdFreetextHelp();
    std::string execute(std::vector<std::string> arguments);

private:
    void addAtomByCubeCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg);
    bool addAtomBySphereCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg);
    C3DVector calcDirVecFromBond(const CAtom* at1, const CAtom* at2, double angleAroundBond, double angleDirBond, double len, int frame) const;
    C3DVector calcDirVecFromAngle(const CAtom* at1, const CAtom* at2, const CAtom* at3,
                                  double angleAroundBond, double angleDirBond, double len, bool& dihedralNotDefined, int frame) const;

private:
    std::string lastError_;
};
