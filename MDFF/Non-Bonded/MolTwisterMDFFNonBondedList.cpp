#include "MolTwisterMDFFNonBondedList.h"

CMDFFNonBondedList::CMDFFNonBondedList()
{
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJ>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJ1208>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_Buck>());
    registeredForceFieldTypes_.emplace_back(std::make_shared<CMDFFNonBonded_LJBuck>());
}

CMDFFNonBondedList::~CMDFFNonBondedList()
{
}

void CMDFFNonBondedList::add(const CMDFFNonBonded& nonBonded)
{
    std::shared_ptr<std::vector<int>> indices = indexFromNames(nonBonded.getAtomInBond(0), nonBonded.getAtomInBond(1));
    appendToList(nonBonded_, *indices, nonBonded);
}

std::shared_ptr<std::vector<int>> CMDFFNonBondedList::indexFromNames(std::string atom1, std::string atom2)
{
    auto indices = std::make_shared<std::vector<int>>();

    for(int i=0; i<(int)nonBonded_.size(); i++)
    {
        if(atom1 == nonBonded_[i]->getAtomInBond(0))
        {
            if(atom2 == nonBonded_[i]->getAtomInBond(1))
            {
                indices->emplace_back(i);
                continue;
            }
        }

        if(atom1 == nonBonded_[i]->getAtomInBond(1))
        {
            if(atom2 == nonBonded_[i]->getAtomInBond(0))
            {
                indices->emplace_back(i);
                continue;
            }
        }
    }

    return indices;
}
