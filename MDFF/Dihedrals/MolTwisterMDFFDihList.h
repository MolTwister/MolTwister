#pragma once
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFDih.h"
#include "MolTwisterMDFFDih_Fourier4t.h"
#include "MolTwisterMDFFDih_Harm.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFDihList : public CMDFFList<CMDFFDih>
{
public:
    CMDFFDihList();
    ~CMDFFDihList();
    
public:
    void add(const CMDFFDih& dih);
    void del(int index) { dihedrals_.erase(dihedrals_.begin() + index); }
    CMDFFDih* get(int index) const { return dihedrals_[index].get(); }
    int size() const { return (int)dihedrals_.size(); }
    void empty() { dihedrals_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string atom1, std::string atom2, std::string atom3, std::string atom4) const;

    int getNumRegisteredFFTypes() { return (int)registeredForceFieldTypes_.size(); }
    CMDFFDih* getRegisteredFFType(int index) { return registeredForceFieldTypes_[index].get(); }
    
private:
    std::vector<std::shared_ptr<CMDFFDih>> dihedrals_;
    std::vector<std::shared_ptr<CMDFFDih>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
