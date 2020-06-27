#pragma once
#include <stdio.h>
#include "MDFF/MolTwisterMDFFList.h"
#include "MolTwisterMDFFAngle.h"
#include "MolTwisterMDFFAngle_Harm.h"
#include "MolTwisterMDFFAngle_Class2.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFAngleList : public CMDFFList<CMDFFAngle>
{
public:
    CMDFFAngleList();
    ~CMDFFAngleList();
    
public:
    void add(const CMDFFAngle& angle);
    void del(int index) { angles_.erase(angles_.begin() + index); }
    CMDFFAngle* get(int index) const { return angles_[index].get(); }
    int size() const { return (int)angles_.size(); }
    void empty() { angles_.clear(); }
    std::shared_ptr<std::vector<int>> indexFromNames(std::string aAtom1, std::string atom2, std::string atom3) const;
    
    int getNumRegisteredFFTypes() const { return (int)registeredForceFieldTypes_.size(); }
    CMDFFAngle* getRegisteredFFType(int index) const { return registeredForceFieldTypes_[index].get(); }

private:
    std::vector<std::shared_ptr<CMDFFAngle>> angles_;
    std::vector<std::shared_ptr<CMDFFAngle>> registeredForceFieldTypes_;
};

END_CUDA_COMPATIBLE()
