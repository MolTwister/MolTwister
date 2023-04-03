//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#pragma once
#include <vector>
#include <map>
#include <string>
#include <memory>
#include "Utilities/Serializer.h"
#include "Utilities/3DVector.h"
#include "Utilities/3DRect.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CAtom
{
public:
    enum EPBCDir { dirX=0, dirY=1, dirZ=2 };

public:
    class CLabel
    {
    public:
        CLabel() = default;
        CLabel(const std::string& name, const C3DVector& displacement, const C3DVector& color) { name_ = name; displacement_ = displacement; color_ = color; }

    public:
        std::string getName() const { return name_; }
        C3DVector getDisplacement() const { return displacement_; }
        C3DVector getColor() const { return color_; }
        void serialize(CSerializer& io, bool saveToStream);

    private:
        std::string name_;
        C3DVector displacement_;
        C3DVector color_;
    };

private:
    class C1to4Conn
    {
    public:
        C1to4Conn() { atom_ = nullptr; numBondsAway_ = -1; }
        
    public:
        CAtom* atom_;
        char numBondsAway_;
    };

    class CBondDest
    {
    public:
        CBondDest() = default;
        CBondDest(CAtom* atom) { atom_ = atom; }

    public:
        void setDest(CAtom* atom) { atom_ = atom; }
        CAtom* getDest() const { return atom_; }
        void setAsDoubleBond(const bool& doubleBond) { isDoubleBond_ = doubleBond; }
        bool isDoubleBond() const { return isDoubleBond_; }

    private:
        CAtom* atom_ = nullptr;
        bool isDoubleBond_ = false;
    };

public:
    CAtom();
    CAtom(double X, double Y, double Z, std::string ID, int atomIndex);
    CAtom(const C3DVector& R, std::string ID, int atomIndex);
    CAtom(const CAtom& src);
    virtual ~CAtom();

public:
    void serialize(CSerializer& io, bool saveToStream, const std::vector<std::shared_ptr<CAtom>>* newAtomList=nullptr);
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
    CAtom* getBondDest(int index) const { return bonds_[index].getDest(); }
    void setAsdoubleBond(int index, bool doubleBond) { bonds_[index].setAsDoubleBond(doubleBond); }
    bool isDoubleBond(int index) const { return bonds_[index].isDoubleBond(); }
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
    void setAtomLabel(const CLabel& label) { atomLabel_ = label; }
    void setAtomLabel(const std::string& name, const C3DVector& displacement, const C3DVector& color) { setAtomLabel(CLabel(name, displacement, color)); }
    CLabel getAtomLabel() const { return atomLabel_; }
    void addBondLabel(CAtom* bondDest, const CLabel& label) { bondLabels_[bondDest] = label; }
    void addBondLabel(CAtom* bondDest, const std::string& name, const C3DVector& displacement, const C3DVector& color) { addBondLabel(bondDest, CLabel(name, displacement, color)); }
    CLabel getBondLabel(CAtom* bondDest) const;

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
    double Q_ = 0.0;
    double sigma_ = 0.0;
    double m_ = 0.0;
    std::string resname_;
    bool isMobile_ = true;
    bool ignoreBondFrom_ = false;

private:
    std::string ID_;
    std::vector<CBondDest> bonds_;
    CLabel atomLabel_;
    std::map<CAtom*, CLabel> bondLabels_;
    int atomIndex_ = -1;
    int molIndex_ = -1;
    bool isSelected_ = false;
    std::vector<C1to4Conn> listOf1to4Connections_;
};

END_CUDA_COMPATIBLE()
