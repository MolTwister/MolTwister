#pragma once
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "Utilities/3DVector.h"
#include "Utilities/3DRect.h"

class CAtom
{
public:
    enum EPBCDir { dirX=0, dirY=1, dirZ=2 };
    
    
private:
    class C1to4Conn
    {
    public:
        C1to4Conn() { atom_ = nullptr; numBondsAway_ = -1; }
        
    public:
        CAtom* atom_;
        char numBondsAway_;
    };
    
    
public:
    CAtom();
    CAtom(double X, double Y, double Z, std::string ID, int atomIndex);
    CAtom(const C3DVector& R, std::string ID, int atomIndex);
    CAtom(const CAtom& src);
    virtual ~CAtom();
    
    
public:
    int addFrame();
    int deleteFrame(int frame);
    void setPos(double X, double Y, double Z, int frame) { r_[frame].set(X, Y, Z); }
    void setPos(const C3DVector& R, int frame) { r_[frame] = R;  }
    void setID(std::string ID) { ID_ = ID; }
    std::string getID() const { return ID_; }
    void setAtomIndex(int index) { atomIndex_ = index; }
    int getAtomIndex() const { return atomIndex_; }
    void setMolIndex(int index) { molIndex_ = index; }
    int getMolIndex() const { return molIndex_; }
    int searchForNonNegMolIndexInAttachedAtoms() const;
    std::shared_ptr<std::vector<CAtom*>> searchForLeafAtomsConnToBond(int bondDest) const;
    int attachBondTo(CAtom* atom);
    int deleteBondTo(CAtom* atom);
    void deleteAllBonds();
    int getNumBonds() const { return (int)bonds_.size(); }
    int getBondDestIndex(const CAtom* atom) const;
    CAtom* getBondDest(int index) const { return bonds_[index]; }
    double getDistanceTo(const CAtom* atom, int frame) const;
    double getDistanceTo2(const CAtom* atom, int frame) const;
    double getDistanceTo2UsingPBC(const CAtom* atom, int frame, C3DRect pbc) const;
    bool isSelected() const { return isSelected_; }
    void select(bool select=true) { isSelected_ = select; }
    void buildListOf1to4Connections();
    int getBondSepTo(const CAtom* atom) const;
    int detectLongestBond(int frame, double& lenghtX, double& lenghtY, double& lenghtZ) const;
    void findAtomsInMolecule(std::vector<CAtom*>& atomsAtPBCBdry1, std::vector<CAtom*>& atomsAtPBCBdry2, C3DRect pbc, EPBCDir pbcDir, int frame);
    void findAtomsInMolecule(std::vector<CAtom*>& atomsInMolecule, int frame);
    int getNum1to4Connections() const { return (int)listOf1to4Connections_.size(); }
    CAtom* get1to4Connection(int index, char& numBondsAway) const { numBondsAway = listOf1to4Connections_[index].numBondsAway_; return listOf1to4Connections_[index].atom_; }
    void copyIntrinsicAtomProperties(const CAtom& src);
    
private:
    void searchForLeafAtomsConnToBond(int bondDest, std::vector<CAtom*>& leafAtoms) const;
    int searchForNonNegMolIndexInAttachedAtoms(std::vector<const CAtom*>& visitedAtoms) const;
    void searchForLeafAtomsConnToBond(const CAtom* bondAt1, const CAtom* bondAt2, std::vector<CAtom*>& leafAtoms) const;
    void copy(const CAtom& src);
    void switchBdry(int& currBdry) const { if(currBdry == 0) currBdry = 1; else currBdry = 0; }
    void findAtomsInMolecule(std::vector<CAtom*>* atomsAtPBCBdry[], std::map<CAtom*,int>& visitedAtoms, int currBdry, double distFromCaller, C3DRect& pbc, EPBCDir pbcDir, int frame);
    void findAtomsInMolecule(std::vector<CAtom*>* atomsInMolecule, std::map<CAtom*,int>& visitedAtoms, int frame);

    
public:
    std::vector<C3DVector> r_;
    double Q_;
    double sigma_;
    double m_;
    std::string resname_;
    bool isMobile_;
    bool ignoreBondFrom_;
    
private:
    std::string ID_;
    std::vector<CAtom*> bonds_;
    int atomIndex_;
    int molIndex_;
    bool isSelected_;
    std::vector<C1to4Conn> listOf1to4Connections_;
};
