#include <iostream>
#include <float.h>
#include <climits>
#include "MolTwisterState.h"


CMolTwisterState::CMolTwisterState()
{
    view3D_ = nullptr;
    currentFrame_ = 0;
}

CMolTwisterState::~CMolTwisterState()
{
    purgeAtomsList();
    purgeGLObjectList();
    purgeVariableList();
}

int CMolTwisterState::deleteFrame(int frame)
{
    int newNumFrames;
    int numFrames = INT_MAX;
    
    if(atoms_.size() == 0) numFrames = 0;
    
    for(int i=0; i<atoms_.size(); i++)
    {
        newNumFrames = atoms_[i]->deleteFrame(frame);
        if(currentFrame_ >= newNumFrames) currentFrame_ = newNumFrames-1;
        if(newNumFrames < numFrames) numFrames = newNumFrames;
    }
    
    return numFrames;
}

void CMolTwisterState::purgeAtomsList()
{
    atoms_.clear();
}

void CMolTwisterState::purgeGLObjectList()
{    
    glObjects_.clear();
}

void CMolTwisterState::purgeVariableList()
{    
    variables_.clear();
}

void CMolTwisterState::purgeFrames(bool keepFirstFrame)
{    
    for(int i=0; i<atoms_.size(); i++)
    {
        C3DVector vectorKeep;

        if(keepFirstFrame && (atoms_[i]->r_.size() > 0))
            vectorKeep = atoms_[i]->r_[0];
        
        atoms_[i]->r_.clear();
        
        if(keepFirstFrame)
            atoms_[i]->r_.emplace_back(vectorKeep);
    }
}

int CMolTwisterState::addFrame()
{
    int newIndex = -2, prev=0;
    bool indicesNotMatching = false;
    
    for(int i=0; i<atoms_.size(); i++)
    {
        newIndex = atoms_[i]->addFrame();
        if(i > 0)
        {
            if(newIndex != prev) indicesNotMatching = true;
        }
        
        prev = newIndex;
    }
    
    if(indicesNotMatching) return -1;
    
    return newIndex;
}

int CMolTwisterState::addGLObject(CGLObject& glObject)
{
    if(glObject.getType() == CGLObject::objLine)
    {
        glObjects_.emplace_back(std::make_shared<CGLObjectLine>(glObject));
    }
    
    return (int)glObjects_.size()-1;
}

void CMolTwisterState::padFrames(int atomIndex, double X, double Y, double Z)
{
    int atInd = atomIndex - 1;
    if(atInd >= 0)
    {
        for(int i=1; i<(int)atoms_[atInd]->r_.size(); i++)
        {
            atoms_[atomIndex]->r_.emplace_back(C3DVector(X, Y, Z));
        }
    }
}

int CMolTwisterState::addAtom(double X, double Y, double Z, std::string ID)
{
    atoms_.emplace_back(std::make_shared<CAtom>(X, Y, Z, ID, (int)atoms_.size()));
    padFrames((int)atoms_.size()-1, X, Y, Z);
    
    return (int)atoms_.size()-1;
}

int CMolTwisterState::addAtom(CAtom& atom)
{
    atom.setAtomIndex((int)atoms_.size());
    atoms_.emplace_back(std::make_shared<CAtom>(atom));
    if(atom.r_.size() > 0)
        padFrames((int)atoms_.size()-1, atom.r_[0].x_, atom.r_[0].y_, atom.r_[0].z_);
    else
        padFrames((int)atoms_.size()-1, 0.0, 0.0, 0.0);
    
    return (int)atoms_.size()-1;
}

void CMolTwisterState::setAtomCoordinates(int frame, int atom, double X, double Y, double Z)
{
    if((atom < atoms_.size()) && (frame < atoms_[atom]->r_.size()))
    {
        atoms_[atom]->r_[frame] = C3DVector(X, Y, Z);
    }
}

int CMolTwisterState::deleteAtom(int index)
{
    CAtom* atomToDelete;
    
    if(atoms_.size() <= 0) return 0;
    if(index >= atoms_.size()) return (int)atoms_.size();
    
    atomToDelete = atoms_[index].get();
    if(atomToDelete)
    {
        for(int i=0; i<atomToDelete->getNumBonds(); i++)
        {
            CAtom* atomBondedToAtomToDel = atomToDelete->getBondDest(i);
            if(atomBondedToAtomToDel)
                atomBondedToAtomToDel->deleteBondTo(atomToDelete);
        }
    }
    
    atoms_.erase(atoms_.begin() + index);
    reassignAtomIndices();
    
    return (int)atoms_.size();
}

void CMolTwisterState::reassignAtomIndices()
{
    for(int i=0; i<atoms_.size(); i++)
    {
        atoms_[i]->setAtomIndex(i);
    }
}

void CMolTwisterState::genMolIndices()
{
    int molIndex;
    int nextMolIndex = 0;
    
    for(int i=0; i<atoms_.size(); i++)
    {
        atoms_[i]->setMolIndex(-1);
    }
    
    for(int i=0; i<atoms_.size(); i++)
    {
        molIndex = atoms_[i]->searchForNonNegMolIndexInAttachedAtoms();
        
        if(molIndex < 0) atoms_[i]->setMolIndex(nextMolIndex++);
        else             atoms_[i]->setMolIndex(molIndex);
    }
}

int CMolTwisterState::addVariable(CVar& variable)
{
    CVar* existingVar;
    std::string text;
    int variableIndex;
    
    if(variable.getType() == CVar::typeEmpty) return (int)variables_.size();
    
    variable.getName(text);
    existingVar = getVariable(text.data(), variableIndex);
    if(existingVar)
    {
        if(variable.getType() == CVar::typeAtom)
        {
            variables_[variableIndex] = std::make_shared<CVarAtom>(variable);
        }
        else if(variable.getType() == CVar::typeBond)
        {
            variables_[variableIndex] = std::make_shared<CVarBond>(variable);
        }
        else if(variable.getType() == CVar::typeAngle)
        {
            variables_[variableIndex] = std::make_shared<CVarAngle>(variable);
        }
        else if(variable.getType() == CVar::typeDihedral)
        {
            variables_[variableIndex] = std::make_shared<CVarDihedral>(variable);
        }
        else
        {
            printf("Error: Attempted to modify to unknown variable type!");
        }
        
        return variableIndex;
    }
    
    else
    {        
        if(variable.getType() == CVar::typeAtom)
        {
            variables_.emplace_back(std::make_shared<CVarAtom>(variable));
        }
        else if(variable.getType() == CVar::typeBond)
        {
            variables_.emplace_back(std::make_shared<CVarBond>(variable));
        }
        else if(variable.getType() == CVar::typeAngle)
        {
            variables_.emplace_back(std::make_shared<CVarAngle>(variable));
        }
        else if(variable.getType() == CVar::typeDihedral)
        {
            variables_.emplace_back(std::make_shared<CVarDihedral>(variable));
        }
        else
        {
            printf("Error: Attempted to add unknown variable type!");
        }
        
        return (int)variables_.size()-1;
    }
}

CVar* CMolTwisterState::getVariable(std::string name, int& variableIndex)
{
    std::string text;
    
    for(int i=0; i<variables_.size(); i++)
    {
        variables_[i]->getName(text);
        if(text == name)
        {
            variableIndex = i;
            return variables_[i].get();
        }
    }
    
    variableIndex = -1;
    
    return nullptr;
}

bool CMolTwisterState::saveCoordinates(int frame)
{
    if(frame < 0) return false;
    
    savedCoordinates_.clear();
    for(int i=0; i<atoms_.size(); i++)
    {
        if(frame >= atoms_[i]->r_.size()) return false;
        savedCoordinates_.emplace_back(atoms_[i]->r_[frame]);
    }
    
    return true;
}

void CMolTwisterState::retrieveSavedCoordinates(int frame)
{
    for(int i=0; i<savedCoordinates_.size(); i++)
    {
        if(i < atoms_.size())
        {
            atoms_[i]->r_[frame] = savedCoordinates_[i];
        }
    }
}

void CMolTwisterState::searchForAtomTypes(std::vector<std::string>& atomTypes, std::vector<std::string>* resnames)
{
    bool alreadyInList;
    
    atomTypes.clear();
    
    for(int i=0; i<atoms_.size(); i++)
    {
        // Get ID (i.e. atom name)
        std::string ID = atoms_[i]->getID();
        
        // Have we concidered this ID before?
        alreadyInList = false;
        for(int j=0; j<atomTypes.size(); j++)
        {
            if(ID == atomTypes[j]) alreadyInList = true;
        }
        
        // If not alredy in list, then add it
        if(!alreadyInList)
        {
            atomTypes.emplace_back(ID);
            if(resnames) resnames->emplace_back(atoms_[i]->resname_);
        }
    }
}

void CMolTwisterState::getMapAtomPtrToIndex(std::map<CAtom*,int>& map) const
{
    map.clear();
    
    for(int i=0; i<atoms_.size(); i++)
    {
        map[atoms_[i].get()] = i;
    }
}

int CMolTwisterState::atomTypeToTypeIndex(const std::vector<std::string>& atomTypes, std::string type)
{
    int atomTypeIndex = -1;
    
    for(int i=0; i<atomTypes.size(); i++)
    {
        if(type == atomTypes[i])
        {
            atomTypeIndex = i;
            break;
        }
    }
    
    return atomTypeIndex;
}

bool CMolTwisterState::getTypeIndexAtomArray(const std::vector<std::string>& atomTypes, std::vector<int>& atomTypeIndices) const
{
    for(int i=0; i<atoms_.size(); i++)
    {
        bool found = false;
        std::string ID = atoms_[i]->getID();
        for(int j=0; j<atomTypes.size(); j++)
        {
            if(ID == atomTypes[j])
            {
                atomTypeIndices.emplace_back(j);
                found = true;
                break;
            }
        }
        
        if(!found)
        {
            printf("Error: could not locate atom type %s in list of atom types!", ID.data());
            return false;
        }
    }
    
    return true;
}

CAtom* CMolTwisterState::getFirstOccurenceOf(std::string ID) const
{
    CAtom* atomPtr = nullptr;

    for(int i=0; i<atoms_.size(); i++)
    {
        std::string cmpID = atoms_[i]->getID();
        if(ID == cmpID)
        {
            atomPtr = atoms_[i].get();
            break;
        }
    }
    
    return atomPtr;
}

C3DRect CMolTwisterState::calcBoundingBox(int frame, const std::vector<int>& atomIndices) const
{
    C3DRect boundBox;
    
    boundBox.rLow_ = C3DVector((double)FLT_MAX, (double)FLT_MAX, (double)FLT_MAX);
    boundBox.rHigh_ = C3DVector(-(double)FLT_MAX, -(double)FLT_MAX, -(double)FLT_MAX);
    
    C3DVector r;
    for(int n=0; n<atomIndices.size(); n++)
    {
        int i = atomIndices[n];
        
        if((i < 0) || (i >= atoms_.size())) continue;
        if(frame >= atoms_[i]->r_.size()) continue;
        
        r = atoms_[i]->r_[frame];
        
        if(r.x_ < boundBox.rLow_.x_) boundBox.rLow_.x_ = r.x_;
        if(r.y_ < boundBox.rLow_.y_) boundBox.rLow_.y_ = r.y_;
        if(r.z_ < boundBox.rLow_.z_) boundBox.rLow_.z_ = r.z_;
        
        if(r.x_ > boundBox.rHigh_.x_) boundBox.rHigh_.x_ = r.x_;
        if(r.y_ > boundBox.rHigh_.y_) boundBox.rHigh_.y_ = r.y_;
        if(r.z_ > boundBox.rHigh_.z_) boundBox.rHigh_.z_ = r.z_;
    }
    
    return boundBox;
}

void CMolTwisterState::getAtomsWithID(std::string ID, std::vector<int>& atomIndices)
{
    atomIndices.clear();
    
    for(int i=0; i<atoms_.size(); i++)
    {
        std::string cmpID = atoms_[i]->getID();
        if(cmpID == ID)
        {
            atomIndices.emplace_back(i);
        }
    }
}

void CMolTwisterState::getAtomsWithID(std::string ID, std::vector<CAtom*>& atomIndices)
{
    atomIndices.clear();
    
    for(int i=0; i<atoms_.size(); i++)
    {
        std::string cmpID = atoms_[i]->getID();
        if(cmpID == ID)
        {
            atomIndices.emplace_back(atoms_[i].get());
        }
    }
}

void CMolTwisterState::getAtomsWithResname(std::string resname, std::vector<int>& atomIndices)
{
    atomIndices.clear();
    
    for(int i=0; i<atoms_.size(); i++)
    {
        if(atoms_[i]->resname_ == resname)
        {
            atomIndices.emplace_back(i);
        }
    }
}

int CMolTwisterState::getAtomIndex(const CAtom* atom) const
{
    for(int i=0; i<atoms_.size(); i++)
    {
        if(atoms_[i].get() == atom) return i;
    }
    
    return -1;
}
