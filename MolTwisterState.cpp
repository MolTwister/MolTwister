//
// Copyright (C) 2024 Richard Olsen.
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

#include <iostream>
#include <float.h>
#include <climits>
#include "Utilities/ASCIIUtility.h"
#include "MolTwisterState.h"

BEGIN_CUDA_COMPATIBLE()

CMolTwisterState::CMolTwisterState()
{
    view3D_ = nullptr;
    currentFrame_ = 0;
    deleteView_ = false;

    // Register variable types
    registeredVariableTypes_.emplace_back(std::make_shared<CVarAtom>());
    registeredVariableTypes_.emplace_back(std::make_shared<CVarBond>());
    registeredVariableTypes_.emplace_back(std::make_shared<CVarAngle>());
    registeredVariableTypes_.emplace_back(std::make_shared<CVarDihedral>());

    // Register GL object types
    registeredGLObjectTypes_.emplace_back(std::make_shared<CGLObjectLine>());
}

CMolTwisterState::~CMolTwisterState()
{
    purgeAtomsList();
    purgeGLObjectList();
    purgeVariableList();

    if(deleteView_) delete view3D_;
}

void CMolTwisterState::serialize(CSerializer& io, bool saveToStream)
{
    // Note: registeredVariableTypes_ is filled in at construction and is used to fulfill
    // serialization and should therefore not itself be seriealized.
    if(saveToStream)
    {
        io << variables_.size();
        int type;
        for(std::shared_ptr<CVar> var : variables_)
        {
            type = (int)var->getType();
            io << type;
            var->serialize(io, saveToStream);
        }

        io << shortcutDirs_.size();
        for(std::string& str : shortcutDirs_)
        {
            io << str;
        }

        io << atoms_.size();
        for(std::shared_ptr<CAtom> atom : atoms_)
        {
            atom->serialize(io, saveToStream);
        }

        io << glObjects_.size();
        for(std::shared_ptr<CGLObject> glObj : glObjects_)
        {
            type = (int)glObj->getType();
            io >> type;
            glObj->serialize(io, saveToStream);
        }

        io << savedCoordinates_.size();
        for(C3DVector vec : savedCoordinates_)
        {
            vec.serialize(io, saveToStream);
        }

        mdFFNonBondedList_.serialize(io, saveToStream);
        mdFFBondList_.serialize(io, saveToStream);
        mdFFAngleList_.serialize(io, saveToStream);
        mdFFDihList_.serialize(io, saveToStream);
        defaultAtProp_.serialize(io, saveToStream);
        view3D_->serialize(io, saveToStream);

        io << currentFrame_;
    }
    else
    {
        int type;
        size_t size;

        io >> size;
        variables_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            io >> type;
            for(std::shared_ptr<CVar>& item : registeredVariableTypes_)
            {
                if(int(item->getType()) == type)
                {
                    variables_[i] = item->createCopy();
                    variables_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }

        io >> size;
        shortcutDirs_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string str;
            io >> str;
            shortcutDirs_[i] = str;
        }

        io >> size;
        atoms_.resize(size);
        for(size_t i=0; i<size; i++) atoms_[i] = std::make_shared<CAtom>();
        for(size_t i=0; i<size; i++)
        {
            atoms_[i]->serialize(io, saveToStream, &atoms_);
        }

        io << glObjects_.size();
        io >> size;
        glObjects_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            io >> type;
            for(std::shared_ptr<CGLObject>& item : registeredGLObjectTypes_)
            {
                if(int(item->getType()) == type)
                {
                    glObjects_[i] = item->createCopy();
                    glObjects_[i]->serialize(io, saveToStream);
                    break;
                }
            }
        }

        io >> size;
        savedCoordinates_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            C3DVector vec;
            vec.serialize(io, saveToStream);
            savedCoordinates_[i] = vec;
        }

        mdFFNonBondedList_.serialize(io, saveToStream);
        mdFFBondList_.serialize(io, saveToStream);
        mdFFAngleList_.serialize(io, saveToStream);
        mdFFDihList_.serialize(io, saveToStream);
        defaultAtProp_.serialize(io, saveToStream);

        deleteView_ = true;
        view3D_ = new C3DView(0, nullptr);
        view3D_->serialize(io, saveToStream, &atoms_, &glObjects_, &currentFrame_, &defaultAtProp_);

        io >> currentFrame_;
    }
}

int CMolTwisterState::deleteFrame(int frame)
{
    int newNumFrames;
    int numFrames = INT_MAX;
    
    if(atoms_.size() == 0) numFrames = 0;
    
    for(int i=0; i<(int)atoms_.size(); i++)
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
    for(int i=0; i<(int)atoms_.size(); i++)
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
    
    for(int i=0; i<(int)atoms_.size(); i++)
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
    bool foundGLType = false;
    for(std::shared_ptr<CGLObject>& item : registeredGLObjectTypes_)
    {
        if(glObject.getType() == item->getType())
        {
            glObjects_.emplace_back(item->createCopy(glObject));
            foundGLType = true;
            break;
        }
    }

    if(!foundGLType)
    {
        printf("Error: Attempted to add unknown GL object type!");
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

std::shared_ptr<std::vector<std::pair<size_t, int>>> CMolTwisterState::getDeletedAtomsChunkSinceLastChunkTransfer() const
{
    auto deletedElements = std::make_shared<std::vector<std::pair<size_t, int>>>();

    // Detect deleted elements by comparing the list already retrieved by the original / modified atom list
    // Note! We are here assuming that the atoms will keep their same order in both lists
    size_t j = 0;
    size_t sizeBeforeChange = atomsRetrievedAsChunks_.size();
    for(size_t i=0; i<sizeBeforeChange; i++)
    {
        std::shared_ptr<CAtom> atomFromBeforeChange = atomsRetrievedAsChunks_[i];
        std::shared_ptr<CAtom> atomFromAfterChange = nullptr;
        if(j < sizeBeforeChange)
        {
            atomFromAfterChange = atoms_[j];
        }

        if(atomFromBeforeChange == atomFromAfterChange)
        {
            j++;
        }
        else
        {
            deletedElements->emplace_back(std::pair<size_t, int>(i, atomFromBeforeChange->getUniqueID()));
        }
    }

    return deletedElements;
}

std::shared_ptr<std::vector<std::shared_ptr<CAtom>>> CMolTwisterState::getAddedAtomsChunkSinceLastChunkTransfer() const
{
    auto addedElements = std::make_shared<std::vector<std::shared_ptr<CAtom>>>();

    // Detect added elements by comparing the list already retrieved by the original / modified atom list
    // Note! We are here assuming that the atoms will keep their same order in both lists
    size_t j = 0;
    size_t sizeAfterChange = atoms_.size();
    size_t sizeBeforeChange = atomsRetrievedAsChunks_.size(); // :TODO: Should we rename to atomsTransferredAsChunks_?
    for(size_t i=0; i<sizeAfterChange; i++)
    {
        std::shared_ptr<CAtom> atomFromAfterChange = atoms_[i];
        if(j >= sizeBeforeChange)
        {
            addedElements->emplace_back(atomFromAfterChange);
        }
        else
        {
            std::shared_ptr<CAtom> atomFromBeforeChange = atomsRetrievedAsChunks_[j++];
            while(atomFromBeforeChange && (atomFromBeforeChange != atomFromAfterChange))
            {
                atomFromBeforeChange = (j < sizeBeforeChange) ? atomsRetrievedAsChunks_[j] : nullptr;
                j++;
            }
        }
    }

    return addedElements;
}

template<class T>
int CMolTwisterState::getChunkCountOfArray(const std::vector<T>& array, const int& numAtomsInChunk) const
{
    return array.size() / numAtomsInChunk + ((array.size() % numAtomsInChunk) ? 1 : 0);
}

template<class T>
void CMolTwisterState::getChunkOfArray(const std::vector<T>& array, const int& numItemsInChunk, const int& chunkIndex,
                                     std::function<void(const T item, const int& index, const bool& firstItemInChunk)> onItemInChunk) const
{
    const int startItemIndex = chunkIndex * numItemsInChunk;
    const int nextItemIndex = startItemIndex + numItemsInChunk;
    const int nextItemIndexAjusted = (nextItemIndex > (int)array.size()) ? (int)array.size() : nextItemIndex;

    for(int i=startItemIndex; i<nextItemIndexAjusted; i++)
    {
        const T item = array[i];
        onItemInChunk(item, i, i == startItemIndex);
    }
}

template<class T>
std::string CMolTwisterState::getChunkOfArrayAsJson(const std::vector<T>& array, const int& numItemsInChunk, const int& chunkIndex, const std::string& jsonObjName,
                                                    std::function<std::string(const T item, const int& index)> onRequestListJsonItem, std::vector<T>* arrayOfRetrievedItemsInChunk) const
{
    std::function<bool(const T, const std::vector<T>&)> alreadyExits = [](const T item, const std::vector<T>& list)
    {
        for(const T& listItem : list)
        {
            if(listItem == item) return true;
        }
        return false;
    };

    std::string jsonItemList = R"({")";
    jsonItemList+= jsonObjName;
    jsonItemList+= R"(":[)";

    std::function<void(const T, const int&, const bool&)> onItemInChunk
        = [&](const T item, const int& index, const bool& firstItemInChunk)
    {
        if(!firstItemInChunk) jsonItemList+= ",";
        jsonItemList+= "{";
        jsonItemList+= onRequestListJsonItem(item, index);
        jsonItemList+= "}";

        if(arrayOfRetrievedItemsInChunk && !alreadyExits(item, *arrayOfRetrievedItemsInChunk)) arrayOfRetrievedItemsInChunk->emplace_back(item);
    };
    getChunkOfArray(array, numItemsInChunk, chunkIndex, onItemInChunk);
    return jsonItemList + "]}";
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
    if((atom < (int)atoms_.size()) && (frame < (int)atoms_[atom]->r_.size()))
    {
        atoms_[atom]->r_[frame] = C3DVector(X, Y, Z);
    }
}

int CMolTwisterState::deleteAtom(int index)
{
    CAtom* atomToDelete;
    
    if(atoms_.size() <= 0) return 0;
    if(index >= (int)atoms_.size()) return (int)atoms_.size();
    
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
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        atoms_[i]->setAtomIndex(i);
    }
}

void CMolTwisterState::genMolIndices()
{
    int molIndex;
    int nextMolIndex = 0;
    
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        atoms_[i]->setMolIndex(-1);
    }
    
    for(int i=0; i<(int)atoms_.size(); i++)
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
        bool foundVariableType = false;
        for(std::shared_ptr<CVar>& item : registeredVariableTypes_)
        {
            if(variable.getType() == item->getType())
            {
                variables_[variableIndex] = item->createCopy(variable);
                foundVariableType = true;
                break;
            }
        }

        if(!foundVariableType)
        {
            printf("Error: Attempted to modify to unknown variable type!");
        }
        
        return variableIndex;
    }
    
    else
    {
        bool foundVariableType = false;
        for(std::shared_ptr<CVar>& item : registeredVariableTypes_)
        {
            if(variable.getType() == item->getType())
            {
                variables_.emplace_back(item->createCopy(variable));
                foundVariableType = true;
                break;
            }
        }

        if(!foundVariableType)
        {
            printf("Error: Attempted to add unknown variable type!");
        }
        
        return (int)variables_.size()-1;
    }
}

CVar* CMolTwisterState::getVariable(std::string name, int& variableIndex)
{
    std::string text;
    
    for(int i=0; i<(int)variables_.size(); i++)
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
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        if(frame >= (int)atoms_[i]->r_.size()) return false;
        savedCoordinates_.emplace_back(atoms_[i]->r_[frame]);
    }
    
    return true;
}

void CMolTwisterState::retrieveSavedCoordinates(int frame)
{
    for(int i=0; i<(int)savedCoordinates_.size(); i++)
    {
        if(i < (int)atoms_.size())
        {
            atoms_[i]->r_[frame] = savedCoordinates_[i];
        }
    }
}

void CMolTwisterState::searchForAtomTypes(std::vector<std::string>& atomTypes, std::vector<std::string>* resnames)
{
    bool alreadyInList;
    
    atomTypes.clear();
    
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        // Get ID (i.e. atom name)
        std::string ID = atoms_[i]->getID();
        
        // Have we concidered this ID before?
        alreadyInList = false;
        for(int j=0; j<(int)atomTypes.size(); j++)
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
    
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        map[atoms_[i].get()] = i;
    }
}

int CMolTwisterState::atomTypeToTypeIndex(const std::vector<std::string>& atomTypes, std::string type)
{
    int atomTypeIndex = -1;
    
    for(int i=0; i<(int)atomTypes.size(); i++)
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
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        bool found = false;
        std::string ID = atoms_[i]->getID();
        for(int j=0; j<(int)atomTypes.size(); j++)
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

    for(int i=0; i<(int)atoms_.size(); i++)
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
    for(int n=0; n<(int)atomIndices.size(); n++)
    {
        int i = atomIndices[n];
        
        if((i < 0) || (i >= (int)atoms_.size())) continue;
        if(frame >= (int)atoms_[i]->r_.size()) continue;
        
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
    
    for(int i=0; i<(int)atoms_.size(); i++)
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
    
    for(int i=0; i<(int)atoms_.size(); i++)
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
    
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        if(atoms_[i]->resname_ == resname)
        {
            atomIndices.emplace_back(i);
        }
    }
}

int CMolTwisterState::getAtomIndex(const CAtom* atom) const
{
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        if(atoms_[i].get() == atom) return i;
    }
    
    return -1;
}

int CMolTwisterState::getAtomsChunkCount(const int& numAtomsInChunk) const
{
    return getChunkCountOfArray<std::shared_ptr<CAtom>>(atoms_, numAtomsInChunk);
}

std::string CMolTwisterState::getAtomsChunkAsJson(const int& numAtomsInChunk, const int& chunkIndex)
{
    std::function<std::string(std::shared_ptr<CAtom>, const int&)> onRequestListJsonItem =
        [&](std::shared_ptr<CAtom> atom, const int& index)
    {
        std::string jsonItem;
        std::tuple<double, double, double> color = CDefaultAtomicProperties().getCPKColor(atom->getID());

        jsonItem+= CASCIIUtility::sprintf(R"("id":"%lu",)", atom->getUniqueID());
        jsonItem+= CASCIIUtility::sprintf(R"("rx":"%g","ry":"%g","rz":"%g",)",
                                               atom->r_[currentFrame_].x_,
                                               atom->r_[currentFrame_].y_,
                                               atom->r_[currentFrame_].z_);
        jsonItem+= CASCIIUtility::sprintf(R"("cr":"%g","cg":"%g","cb":"%g")",
                                               std::get<0>(color),
                                               std::get<1>(color),
                                               std::get<2>(color));
        return jsonItem;
    };

    return getChunkOfArrayAsJson<std::shared_ptr<CAtom>>(atoms_, numAtomsInChunk, chunkIndex, "AtomsChunk", onRequestListJsonItem, &atomsRetrievedAsChunks_);
}

/* :TODO: Remove after testing the above function
//        const CAtom* atom = atoms_[i].get();
std::string CMolTwisterState::getAtomsChunkAsJson(const int& numAtomsInChunk, const int& chunkIndex)
{
    std::string jsonAtomList = R"({"AtomsChunk":[)";
    std::function<void(const CAtom* atom, const int& index, const bool& firstAtomInChunk)> onAtomInChunk
        = [&](const CAtom* atom, const int& index, const bool& firstAtomInChunk)
        {
            if(!firstAtomInChunk) jsonAtomList+= ",";
            jsonAtomList+= "{";
            std::tuple<double, double, double> color = CDefaultAtomicProperties().getCPKColor(atom->getID());

            jsonAtomList+= CASCIIUtility::sprintf(R"("rx":"%g","ry":"%g","rz":"%g",)",
                                                 atom->r_[currentFrame_].x_,
                                                 atom->r_[currentFrame_].y_,
                                                 atom->r_[currentFrame_].z_);

            jsonAtomList+= CASCIIUtility::sprintf(R"("cr":"%g","cg":"%g","cb":"%g")",
                                                 std::get<0>(color),
                                                 std::get<1>(color),
                                                 std::get<2>(color));
            jsonAtomList+= "}";
        };
    getAtomsChunk(numAtomsInChunk, chunkIndex, onAtomInChunk);
    return jsonAtomList + "]}";
}
*/
/*
std::string CMolTwisterState::getAtomsChunkAsJson(const int& numAtomsInChunk, const int& chunkIndex)
{
    const int startAtomIndex = chunkIndex * numAtomsInChunk;
    const int nextAtomIndex = startAtomIndex + numAtomsInChunk;
    const int nextAtomIndexAjusted = (nextAtomIndex > (int)atoms_.size()) ? (int)atoms_.size() : nextAtomIndex;

    std::string jsonAtomList = R"({"AtomsChunk":[)";
    for(int i=startAtomIndex; i<nextAtomIndexAjusted; i++)
    {
        if(i != startAtomIndex) jsonAtomList+= ",";

        jsonAtomList+= "{";
        const CAtom* atom = atoms_[i].get();
        std::tuple<double, double, double> color = CDefaultAtomicProperties().getCPKColor(atom->getID());

        jsonAtomList+= CASCIIUtility::sprintf(R"("rx":"%g","ry":"%g","rz":"%g",)",
                                               atom->r_[currentFrame_].x_,
                                               atom->r_[currentFrame_].y_,
                                               atom->r_[currentFrame_].z_);

        jsonAtomList+= CASCIIUtility::sprintf(R"("cr":"%g","cg":"%g","cb":"%g")",
                                               std::get<0>(color),
                                               std::get<1>(color),
                                               std::get<2>(color));
        jsonAtomList+= "}";
    }

    return jsonAtomList + "]}";
}
*/

void CMolTwisterState::resetRetrievedChunksStatus()
{
    atomsRetrievedAsChunks_.erase(atomsRetrievedAsChunks_.begin(), atomsRetrievedAsChunks_.end());
}

int CMolTwisterState::getDeletedAtomsChunkCount(const int& numAtomsInChunk) const
{
    std::shared_ptr<std::vector<std::pair<size_t, int>>> deletedElements = getDeletedAtomsChunkSinceLastChunkTransfer();
    return getChunkCountOfArray<std::pair<size_t, int>>(*deletedElements, numAtomsInChunk);
}

std::string CMolTwisterState::getDeletedAtomsChunkSinceLastChunkTransferAsJson(const int& numAtomsInChunk, const int& chunkIndex)
{
    if(chunkIndex < 0) return R"({"DeletedIndicesChunk":[]})";

    // Detect deleted elements by comparing the list already retrieved by the original / modified atom list
    // Note! We are here assuming that the atoms will keep their same order in both lists
    std::shared_ptr<std::vector<std::pair<size_t, int>>> deletedElements = getDeletedAtomsChunkSinceLastChunkTransfer();

    // Remove all the detected elements from the state, but start from the highest index to avoid shifting index positions.
    // Note that the highest index is always last in the state list. However, only do this when the last chunk is reveived.
    if(chunkIndex >= (getChunkCountOfArray<std::pair<size_t, int>>(*deletedElements, numAtomsInChunk)-1))
    {
        for(int i=(int)deletedElements->size()-1; i>=0; i--)
        {
            atomsRetrievedAsChunks_.erase(atomsRetrievedAsChunks_.begin() + (*deletedElements)[(size_t)i].first);
        }
    }

    // Create a JSON list of the deleted indices and return it using the defined chunk size.
    // The first element of deletedElementIndex is the index into the array before deleting
    // and the second is the unique id of the atom
    std::function<std::string(std::pair<size_t, int>, const int&)> onRequestListJsonItem =
        [&](std::pair<size_t, int> deletedElement, const int& index)
    {
        const unsigned long uniqueID = deletedElement.second;
        return CASCIIUtility::sprintf(R"("id":"%lu")", uniqueID);
    };

    return getChunkOfArrayAsJson<std::pair<size_t, int>>(*deletedElements, numAtomsInChunk, chunkIndex, "DeletedIndicesChunk", onRequestListJsonItem);
}

int CMolTwisterState::getAddedAtomsChunkCount(const int& numAtomsInChunk) const
{
    std::shared_ptr<std::vector<std::shared_ptr<CAtom>>> addedElements = getAddedAtomsChunkSinceLastChunkTransfer();
    return getChunkCountOfArray<std::shared_ptr<CAtom>>(*addedElements, numAtomsInChunk);
}

std::string CMolTwisterState::getAddedAtomsChunkSinceLastChunkTransferAsJson(const int& numAtomsInChunk, const int& chunkIndex)
{
    if(chunkIndex < 0) return R"({"AddedAtomsChunk":[]})";

    // Detect added elements by comparing the list already retrieved by the original / modified atom list
    // Note! We are here assuming that the atoms will keep their same order in both lists
    std::shared_ptr<std::vector<std::shared_ptr<CAtom>>> addedElements = getAddedAtomsChunkSinceLastChunkTransfer();

    // Add all the added elements to the end of the state list. However, only do this when the last chunk is reveived.
    if(chunkIndex >= (getChunkCountOfArray<std::shared_ptr<CAtom>>(*addedElements, numAtomsInChunk)-1))
    {
        for(const std::shared_ptr<CAtom>& atom : *addedElements)
        {
            atomsRetrievedAsChunks_.emplace_back(atom);
        }
    }

    // Create a JSON list of the added atoms and return it using the defined chunk size
    std::function<std::string(std::shared_ptr<CAtom>, const int&)> onRequestListJsonItem =
        [&](std::shared_ptr<CAtom> atom, const int& index)
    {
        // :TODO: There is a similar chunk of code above, create a common function for this!
        std::string jsonItem;
        std::tuple<double, double, double> color = CDefaultAtomicProperties().getCPKColor(atom->getID());

        jsonItem+= CASCIIUtility::sprintf(R"("id":"%lu",)", atom->getUniqueID());
        jsonItem+= CASCIIUtility::sprintf(R"("rx":"%g","ry":"%g","rz":"%g",)",
                                           atom->r_[currentFrame_].x_,
                                           atom->r_[currentFrame_].y_,
                                           atom->r_[currentFrame_].z_);
        jsonItem+= CASCIIUtility::sprintf(R"("cr":"%g","cg":"%g","cb":"%g")",
                                           std::get<0>(color),
                                           std::get<1>(color),
                                           std::get<2>(color));
        return jsonItem;
    };

    return getChunkOfArrayAsJson<std::shared_ptr<CAtom>>(*addedElements, numAtomsInChunk, chunkIndex, "AddedAtomsChunk", onRequestListJsonItem);
}

END_CUDA_COMPATIBLE()
