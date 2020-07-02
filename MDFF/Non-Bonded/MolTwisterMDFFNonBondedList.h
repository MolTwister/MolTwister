#pragma once
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFNonBonded.h"
#include "MolTwisterMDFFNonBonded_LJ.h"
#include "MolTwisterMDFFNonBonded_LJ1208.h"
#include "MolTwisterMDFFNonBonded_Buck.h"
#include "MolTwisterMDFFNonBonded_LJBuck.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFNonBondedList : public CMDFFList<CMDFFNonBonded>
{
public:
    CMDFFNonBondedList();
    ~CMDFFNonBondedList();
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    void add(const CMDFFNonBonded& nonBonded);
    void del(int index) { nonBonded_.erase(nonBonded_.begin() + index); }
    CMDFFNonBonded* get(int index) const { return nonBonded_[index].get(); }
    int size() const { return (int)nonBonded_.size(); }
    void empty() { nonBonded_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2);
    
    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFNonBonded* getRegisteredFFType(int iIndex) const { return registeredForceFieldTypes_[iIndex].get(); }

private:
    std::vector<std::shared_ptr<CMDFFNonBonded>> nonBonded_;
    std::vector<std::shared_ptr<CMDFFNonBonded>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
