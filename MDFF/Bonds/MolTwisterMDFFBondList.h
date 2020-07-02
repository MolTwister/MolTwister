#pragma once
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFBond.h"
#include "MolTwisterMDFFBond_Harm.h"
#include "MolTwisterMDFFBond_Morse.h"
#include "MolTwisterMDFFBond_LJC.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFBondList : public CMDFFList<CMDFFBond>
{
public:
    CMDFFBondList();
    ~CMDFFBondList();
    
public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    void add(const CMDFFBond& bond);
    void del(int index) { bonds_.erase(bonds_.begin() + index); }
    CMDFFBond* get(int index) const { return bonds_[index].get(); }
    int size() const { return (int)bonds_.size(); }
    void empty() { bonds_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2) const;

    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFBond* getRegisteredFFType(int index) const { return registeredForceFieldTypes_[index].get(); }
    
private:
    std::vector<std::shared_ptr<CMDFFBond>> bonds_;
    std::vector<std::shared_ptr<CMDFFBond>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
