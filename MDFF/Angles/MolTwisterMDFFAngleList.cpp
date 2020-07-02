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

void CMDFFAngleList::serialize(std::stringstream& io, bool saveToStream)
{
    CMDFFList<CMDFFAngle>::serialize(io, saveToStream);

    // Note: we do not serialize registeredForceFieldTypes_, size this is always built in the
    // class constructor. Moreover, it is used to read from the stream
    if(saveToStream)
    {
        io << angles_.size();
        for(std::shared_ptr<CMDFFAngle> item : angles_)
        {
            io << item->getFFType();
            item->serialize(io, saveToStream);
        }
    }
    else
    {
        size_t size;

        io >> size;
        angles_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            std::string ffType;
            io >> ffType;
            for(std::shared_ptr<CMDFFAngle> item : registeredForceFieldTypes_)
            {
                if(item->getFFType() == ffType)
                {
                    std::shared_ptr<CMDFFAngle> copy = item->createCopy();
                    copy->serialize(io, saveToStream);
                    angles_[i] = copy;
                    break;
                }
            }
        }
    }
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
