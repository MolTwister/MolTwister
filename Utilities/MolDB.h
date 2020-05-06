#pragma once
#include <stdio.h>
#include <vector>
#include <string>

class CMolDB
{
public:
    class CAtom
    {
    public:
        CAtom();
        
    public:
        bool parse(std::string line);
        void print() const;
        bool operator==(const CAtom& src) const;

    public:
        std::string name_;
        double x_;
        double y_;
        double z_;
        double mass_;
        double q_;
        double ljSigma_;
        double ljEpsilon_;
        double buckA_;
        double buckRho_;
        double buckC_;
        int atomType_;
        int atomTypeGlobal_;
        bool firstAtomOfThisType_;
    };
    
    class CBond
    {
    public:
        CBond();
        
    public:
        bool parse(std::string line);
        
    public:
        int connection_[2];
        char type_;
        double param1_;
        double param2_;
        double param3_;
        double param4_;
    };

    class CAngle
    {
    public:
        CAngle();
        
    public:
        bool parse(std::string line);
        
    public:
        int connection_[3];
        double param1_;
        double param2_;
    };

    class CDihedral
    {
    public:
        CDihedral();
        
    public:
        bool parse(std::string line, int torsType);
        int getTorsType() const { return torsType_; }
        
    public:
        int connection_[4];
        double param1_;
        double param2_;
        double param3_;
        double param4_;
        double param5_;
        int torsType_;
    };
    
    class CMolecule
    {
    public:
        CMolecule() { numComments_ = 0; flexibleSPC_ = false; }
        
    public:
        bool load(FILE* fileHandle, int offsetAtomTypeGlobal=0);
        int getNumAtoms() const { return (int)atoms_.size(); }
        int getNumDistinctAtomTypes() const;
        int getNumBonds() const { return (int)bonds_.size(); }
        int getNumAngles() const { return (int)angles_.size(); }
        int getNumDihedrals() const { return (int)dihedrals_.size(); }
        CAtom* getAtom(int index) { return &atoms_[index]; }
        const CAtom* getAtom(int index) const { return &atoms_[index]; }
        CBond* getBond(int index) { return &bonds_[index]; }
        const CBond* getBond(int index) const { return &bonds_[index]; }
        CAngle* getAngle(int index) { return &angles_[index]; }
        const CAngle* getAngle(int index) const { return &angles_[index]; }
        CDihedral* getDihedral(int index) { return &dihedrals_[index]; }
        const CDihedral* getDihedral(int index) const { return &dihedrals_[index]; }
        void print() const;
        bool isFlexibleSPC() const { return flexibleSPC_; }
        
    private:
        bool handleTors(FILE* fileHandle, int torsType);
        int searchForIdenticalAtomType(const CAtom& atom, int lastIndex=-1) const;

    private:
        int numComments_;
        std::vector<CAtom> atoms_;
        std::vector<CBond> bonds_;
        std::vector<CAngle> angles_;
        std::vector<CDihedral> dihedrals_;
        bool flexibleSPC_;
    };
    
public:
    CMolDB();
    
public:
    bool addMoleculeType(const char* fileName);
    int getNumMoleculeTypes() const { return (int)moleculeTypes_.size(); }
    CMolecule* getMoleculeType(int molTypeIndex) { return &moleculeTypes_[molTypeIndex]; }
    void purge();
    void print() const;

private:
    std::vector<CMolecule> moleculeTypes_;
};
