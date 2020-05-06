#pragma once
#include <vector>
#include <string>
#include <map>
#include <Python.h>
#include "MolTwisterAtom.h"
#include "3DView/MolTwisterGLObject.h"
#include "MolTwisterVariable.h"
#include "3DView/MolTwister3DView.h"
#include "Utilities/3DRect.h"
#include "DefaultAtomicProperties.h"
#include "MDFF/Non-Bonded/MolTwisterMDFFNonBondedList.h"
#include "MDFF/Bonds/MolTwisterMDFFBondList.h"
#include "MDFF/Angles/MolTwisterMDFFAngleList.h"
#include "MDFF/Dihedrals/MolTwisterMDFFDihList.h"

#define MOLTWISTER_VER "1.4.0"

class CMolTwisterState
{
public:
    CMolTwisterState();
    ~CMolTwisterState();
    
public:
    void purgeAtomsList();
    void purgeGLObjectList();
    void purgeVariableList();
    void purgeFrames(bool keepFirstFrame=true);
    int addFrame();
    int addAtom(double X, double Y, double Z, std::string ID);
    int addAtom(CAtom& atom);
    void setAtomCoordinates(int frame, int atom, double X, double Y, double Z);
    int deleteAtom(int index);
    int deleteFrame(int frame);
    int addGLObject(CGLObject& glObject);
    int addVariable(CVar& variable);
    CVar* getVariable(std::string name, int& variableIndex);
    void reassignAtomIndices();
    void genMolIndices();
    int getCurrFrameIndex() const { return currentFrame_; }
    bool saveCoordinates(int frame);
    void retrieveSavedCoordinates(int frame);
    void searchForAtomTypes(std::vector<std::string>& atomTypes, std::vector<std::string>* resnames=nullptr);
    static int atomTypeToTypeIndex(const std::vector<std::string>& atomTypes, std::string type);
    bool getTypeIndexAtomArray(const std::vector<std::string>& atomTypes, std::vector<int>& atomTypeIndices) const;
    void getMapAtomPtrToIndex(std::map<CAtom*,int>& map) const;
    CAtom* getFirstOccurenceOf(std::string ID) const;
    C3DRect calcBoundingBox(int frame, const std::vector<int>& atomIndices) const;
    void getAtomsWithID(std::string ID, std::vector<int>& atomIndices);
    void getAtomsWithID(std::string ID, std::vector<CAtom*>& atomIndices);
    void getAtomsWithResname(std::string resname, std::vector<int>& atomIndices);
    int getAtomIndex(const CAtom* atom) const;

private:
    void padFrames(int atomIndex, double X, double Y, double Z);
    
public:
    CMDFFNonBondedList mdFFNonBondedList_;
    CMDFFBondList mdFFBondList_;
    CMDFFAngleList mdFFAngleList_;
    CMDFFDihList mdFFDihList_;
    std::vector<std::shared_ptr<CVar>> variables_;
    std::vector<std::string> shortcutDirs_;
    std::vector<std::shared_ptr<CAtom>> atoms_;
    std::vector<std::shared_ptr<CGLObject>> glObjects_;
    std::vector<C3DVector> savedCoordinates_;
    CDefaultAtomicProperties defaultAtProp_;
    C3DView* view3D_;
    int currentFrame_;
};
