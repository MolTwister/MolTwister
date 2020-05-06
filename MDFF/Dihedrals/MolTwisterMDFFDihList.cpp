#include "MolTwisterMDFFDihList.h"

CMDFFDihList::CMDFFDihList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFDih_Fourier4t>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFDih_Harm>());
}

CMDFFDihList::~CMDFFDihList()
{
}

void CMDFFDihList::add(const CMDFFDih& dih)
{    
    std::shared_ptr<std::vector<int>> indices = indexFromNames(dih.getAtomInBond(0), dih.getAtomInBond(1), dih.getAtomInBond(2), dih.getAtomInBond(3));
    appendToList(dihedrals_, *indices, dih);
}

std::shared_ptr<std::vector<int>> CMDFFDihList::indexFromNames(std::string atom1, std::string atom2, std::string atom3, std::string atom4) const
{
    auto indices = std::make_shared<std::vector<int>>();
    
    for(int i=0; i<dihedrals_.size(); i++)
    {
        if(atom1 == dihedrals_[i]->getAtomInBond(0))
        {
            if(atom2 == dihedrals_[i]->getAtomInBond(1))
            {
                if(atom3 == dihedrals_[i]->getAtomInBond(2))
                {
                    if(atom4 == dihedrals_[i]->getAtomInBond(3))
                    {
                        indices->emplace_back(i);
                        continue;
                    }
                }
            }
        }
        
        if(atom1 == dihedrals_[i]->getAtomInBond(3))
        {
            if(atom2 == dihedrals_[i]->getAtomInBond(2))
            {
                if(atom3 == dihedrals_[i]->getAtomInBond(1))
                {
                    if(atom4 == dihedrals_[i]->getAtomInBond(0))
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
