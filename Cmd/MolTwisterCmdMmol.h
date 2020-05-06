#pragma once
#include "MolTwisterCmd.h"
#include "Utilities/MolDB.h"

class CCmdMmol : public CCmd
{
private:
    struct SPairInt
    {
        SPairInt() { }
        SPairInt(std::string s1, std::string s2) { atom1_ = s1; atom2_ = s2; }
        std::string atom1_;
        std::string atom2_;
    };

    struct SAngle
    {
        SAngle() { }
        SAngle(std::string s1, std::string s2, std::string s3) { atom1_ = s1; atom2_ = s2; atom3_ = s3; }
        std::string atom1_;
        std::string atom2_;
        std::string atom3_;
    };

    struct SDih
    {
        SDih() { }
        SDih(std::string s1, std::string s2, std::string s3, std::string s4) { atom1_ = s1; atom2_ = s2; atom3_ = s3; atom4_ = s4; }
        std::string atom1_;
        std::string atom2_;
        std::string atom3_;
        std::string atom4_;
    };
    
public:
    CCmdMmol() = delete;
    CCmdMmol(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "mmol"; }
    std::string getTopLevHelpString() const { return std::string("Apply operations to *.mmol type files"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseToscriptCommand(std::string commandLine, int& arg);
    bool doesDihExist(std::string atom1, std::string atom2, std::string atom3, std::string atom4, const std::vector<SDih>& dihedrals) const;
    bool doesAngleExist(std::string atom1, std::string atom2, std::string atom3, const std::vector<SAngle>& angles) const;
    bool doesIntExist(std::string atom1, std::string atom2, const std::vector<SPairInt>& pairInt) const;
    bool doesAtomExist(std::string atom, const std::vector<std::string>& atoms) const;
    void genVanDerWaals(const CMolDB::CMolecule* molecule1, const CMolDB::CMolecule* molecule2, bool mixGeom) const;
};
