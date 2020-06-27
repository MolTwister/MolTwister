#include "MolTwisterMDFFBondList.h"

BEGIN_CUDA_COMPATIBLE()

CMDFFBondList::CMDFFBondList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_Harm>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_Morse>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFBond_LJC>());
}

CMDFFBondList::~CMDFFBondList()
{
}

void CMDFFBondList::add(const CMDFFBond& bond)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(bond.getAtomInBond(0), bond.getAtomInBond(1));
    appendToList(bonds_, *indices, bond);
}

std::shared_ptr<std::vector<int>> CMDFFBondList::indexFromNames(std::string atom1, std::string atom2) const
{
    auto indices = std::make_shared<std::vector<int>>();
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(atom1 == bonds_[i]->getAtomInBond(0))
        {
            if(atom2 == bonds_[i]->getAtomInBond(1))
            {
                indices->emplace_back(i);
                continue;
            }
        }
        
        if(atom1 == bonds_[i]->getAtomInBond(1))
        {
            if(atom2 == bonds_[i]->getAtomInBond(0))
            {
                indices->emplace_back(i);
                continue;
            }
        }
    }

    return indices;
}

END_CUDA_COMPATIBLE()
