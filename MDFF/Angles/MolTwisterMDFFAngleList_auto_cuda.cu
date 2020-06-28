#include "MolTwisterMDFFAngleList.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFAngleList::CMDFFAngleList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFAngle_Harm>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFAngle_Class2>());
}

CMDFFAngleList::~CMDFFAngleList()
{
}

void CMDFFAngleList::add(const CMDFFAngle& angle)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(angle.getAtomInBond(0), angle.getAtomInBond(1), angle.getAtomInBond(2));
    appendToList(angles_, *indices, angle);
}

std::shared_ptr<std::vector<int>> CMDFFAngleList::indexFromNames(std::string aAtom1, std::string atom2, std::string atom3) const
{
    auto indices = std::make_shared<std::vector<int>>();

    for(int i=0; i<(int)angles_.size(); i++)
    {
        if(aAtom1 == angles_[i]->getAtomInBond(0))
        {
            if(atom2 == angles_[i]->getAtomInBond(1))
            {
                if(atom3 == angles_[i]->getAtomInBond(2))
                {
                    indices->emplace_back(i);
                    continue;
                }
            }
        }
        
        if(angles_[i]->isAngleABCEqualToCBA())
        {
            if(aAtom1 == angles_[i]->getAtomInBond(2))
            {
                if(atom2 == angles_[i]->getAtomInBond(1))
                {
                    if(atom3 == angles_[i]->getAtomInBond(0))
                    {
                        indices->emplace_back(i);
                        continue;
                    }
                }
            }
        }
    }

    return indices;
}

END_CUDA_COMPATIBLE()
