#include <iostream>
#include <math.h>
#include <assert.h>
#include "Utilities/3DRect.h"
#include "MolTwisterAtom.h"

BEGIN_CUDA_COMPATIBLE()

CAtom::CAtom()
{
    atomIndex_ = -1;
    molIndex_ = -1;
    isSelected_ = false;

    r_.emplace_back(C3DVector(0.0, 0.0, 0.0));
    Q_ = 0.0;
    m_ = 0.0;
    sigma_ = 1.0;
    ID_ = "NA";
    isMobile_ = true;
    ignoreBondFrom_ = false;
}

CAtom::CAtom(double X, double Y, double Z, std::string ID, int atomIndex)
{
    atomIndex_ = atomIndex;
    molIndex_ = -1;
    isSelected_ = false;    
    
    r_.emplace_back(C3DVector(0.0, 0.0, 0.0));
    r_[0].set(X, Y, Z);
    Q_ = 0.0;
    m_ = 0.0;
    sigma_ = 1.0;
    ID_ = ID;
    isMobile_ = true;
    ignoreBondFrom_ = false;
}

CAtom::CAtom(const C3DVector &R, std::string ID, int atomIndex)
{
    atomIndex_ = atomIndex;
    molIndex_ = -1;
    isSelected_ = false;
    
    r_.emplace_back(C3DVector(0.0, 0.0, 0.0));
    r_[0] = R;
    Q_ = 0.0;
    m_ = 0.0;
    sigma_ = 1.0;
    ID_ = ID;
    isMobile_ = true;
    ignoreBondFrom_ = false;
}

CAtom::CAtom(const CAtom& src)
{ 
    ID_[0] = '\0';

    copy(src); 
}

CAtom::~CAtom()
{
}

void CAtom::serialize(CSerializer& io, bool saveToStream, const std::vector<std::shared_ptr<CAtom>>* newAtomList)
{
    if(saveToStream)
    {
        io << r_.size();
        for(C3DVector item : r_)
        {
            item.serialize(io, saveToStream);
        }

        io << bonds_.size();
        for(CAtom* atom : bonds_)
        {
            io << atom->atomIndex_;
        }

        io << listOf1to4Connections_.size();
        for(C1to4Conn item : listOf1to4Connections_)
        {
            io << item.atom_->atomIndex_;
            io << item.numBondsAway_;
        }

        io << Q_;
        io << sigma_;
        io << m_;
        io << resname_;
        io << isMobile_;
        io << ignoreBondFrom_;
        io << ID_;
        io << atomIndex_;
        io << molIndex_;
        io << isSelected_;
    }
    else
    {
        assert(newAtomList != nullptr);

        size_t size;

        io >> size;
        r_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            r_[i].serialize(io, saveToStream);
        }

        io >> size;
        bonds_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            int atomIndex;
            io >> atomIndex;
            bonds_[i] = (*newAtomList)[atomIndex].get();
        }

        io >> size;
        listOf1to4Connections_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            C1to4Conn conn;

            int atomIndex;
            io >> atomIndex;
            conn.atom_ = (*newAtomList)[atomIndex].get();
            io >> conn.numBondsAway_;

            listOf1to4Connections_[i] = conn;
        }

        io >> Q_;
        io >> sigma_;
        io >> m_;
        io >> resname_;
        io >> isMobile_;
        io >> ignoreBondFrom_;
        io >> ID_;
        io >> atomIndex_;
        io >> molIndex_;
        io >> isSelected_;
    }
}

int CAtom::addFrame()
{
    if(r_.size() == 0)   r_.emplace_back(C3DVector(0.0, 0.0, 0.0));
    else                r_.emplace_back(r_[r_.size()-1]);
    
    return (int)r_.size() - 1;
}

int CAtom::deleteFrame(int frame)
{
    if(frame < (int)r_.size())
        r_.erase(r_.begin() + frame);

    return (int)r_.size();
}

int CAtom::searchForNonNegMolIndexInAttachedAtoms() const
{
    std::vector<const CAtom*> visitedAtoms;
    int index;
    
    if(molIndex_ > -1) return molIndex_;

    visitedAtoms.emplace_back(this);
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        index = bonds_[i]->searchForNonNegMolIndexInAttachedAtoms(visitedAtoms);
        if(index > -1)
        {
            return index;
        }
    }

    return -1;
}

void CAtom::searchForLeafAtomsConnToBond(int bondDest, std::vector<CAtom*>& leafAtoms) const
{
    CAtom* bondDestPtr = nullptr;
    if((bondDest >= 0) && (bondDest < (int)bonds_.size()))
        bondDestPtr = getBondDest(bondDest);

    leafAtoms.clear();

    bool alreadyVisited;
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        alreadyVisited = false;
        for(int j=0; j<(int)leafAtoms.size(); j++)
        {
            if(bonds_[i] == leafAtoms[j]) alreadyVisited = true;
        }
        if(bonds_[i] == bondDestPtr) alreadyVisited = true;

        if(alreadyVisited) continue;

        leafAtoms.emplace_back(bonds_[i]);
        bonds_[i]->searchForLeafAtomsConnToBond(bondDestPtr, this, leafAtoms);
    }
}

int CAtom::searchForNonNegMolIndexInAttachedAtoms(std::vector<const CAtom*>& visitedAtoms) const
{
    bool alreadyVisited;
    int index;

    if(molIndex_ > -1) return molIndex_;

    visitedAtoms.emplace_back(this);
    if(visitedAtoms.size() > 200)
    {
        return -1;
    }

    for(int i=0; i<(int)bonds_.size(); i++)
    {
        alreadyVisited = false;
        for(int j=0; j<(int)visitedAtoms.size(); j++)
        {
            if(visitedAtoms[j] == bonds_[i]) alreadyVisited = true;
        }
        
        if(alreadyVisited) continue;

        index = bonds_[i]->searchForNonNegMolIndexInAttachedAtoms(visitedAtoms);
        if(index > -1) return index;
    }
    
    return -1;
}

std::shared_ptr<std::vector<CAtom*>> CAtom::searchForLeafAtomsConnToBond(int bondDest) const
{
    auto leafAtoms = std::make_shared<std::vector<CAtom*>>();
    searchForLeafAtomsConnToBond(bondDest, *leafAtoms);
    return leafAtoms;
}

void CAtom::searchForLeafAtomsConnToBond(const CAtom* bondAt1, const CAtom* bondAt2, std::vector<CAtom*>& leafAtoms) const
{
    bool alreadyVisited;
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        alreadyVisited = false;
        for(int j=0; j<(int)leafAtoms.size(); j++)
        {
            if(bonds_[i] == leafAtoms[j]) alreadyVisited = true;
        }
        if(bonds_[i] == bondAt1) alreadyVisited = true;
        if(bonds_[i] == bondAt2) alreadyVisited = true;
        
        if(alreadyVisited) continue;
        
        leafAtoms.emplace_back(bonds_[i]);
        bonds_[i]->searchForLeafAtomsConnToBond(bondAt1, bondAt2, leafAtoms);
    }
}

int CAtom::attachBondTo(CAtom* atom)
{
    bonds_.emplace_back(atom);

    return (int)bonds_.size();
}

int CAtom::deleteBondTo(CAtom* atom)
{
    int index = -1;
    
    if(bonds_.size() < 1) return -1;
    
    // Find element index to erase
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(bonds_[i] == atom)
        {
            index = i;
            break;
        }
    }
    
    // If index was found, then delete
    if(index != -1)
    {
        bonds_.erase(bonds_.begin() + index);
    }

    return (int)bonds_.size();
}

void CAtom::deleteAllBonds()
{
    bonds_.clear();
}

int CAtom::getBondDestIndex(const CAtom* atom) const
{
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(bonds_[i] == atom) return i;
    }
    
    return -1;
}

double CAtom::getDistanceTo(const CAtom* atom, int frame) const
{
    return (r_[frame] - atom->r_[frame]).norm();
}

double CAtom::getDistanceTo2(const CAtom* atom, int frame) const
{
    return (r_[frame] - atom->r_[frame]).norm2();
}

double CAtom::getDistanceTo2UsingPBC(const CAtom* atom, int frame, C3DRect pbc) const
{
    double dx = r_[frame].x_ - atom->r_[frame].x_;
    double dy = r_[frame].y_ - atom->r_[frame].y_;
    double dz = r_[frame].z_ - atom->r_[frame].z_;
    double w;

    dx = (dx < 0.0) ? -dx : dx;
    dy = (dy < 0.0) ? -dy : dy;
    dz = (dz < 0.0) ? -dz : dz;
    
    w = pbc.getWidthX();
    if(dx > (w / 2.0)) dx = w - dx;

    w = pbc.getWidthY();
    if(dy > (w / 2.0)) dy = w - dy;

    w = pbc.getWidthZ();
    if(dz > (w / 2.0)) dz = w - dz;
    
    return (dx*dx + dy*dy + dz*dz);
}

void CAtom::buildListOf1to4Connections()
{
    CAtom* atom1Ptr;
    CAtom* atom2Ptr;
    CAtom* atom3Ptr;
    CAtom* atom4Ptr;
    C1to4Conn connect;
    
    // Go through all atoms that are from 1 to 4 bond lengths
    // from this atom and build a complete list. In each entry
    // make a note of the atomic pointer and how many bond lengths
    // from this atom it is located.
    listOf1to4Connections_.clear();
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        atom1Ptr = bonds_[i];
        if(!atom1Ptr) continue;
        connect.atom_ = atom1Ptr;
        connect.numBondsAway_ = 1;
        listOf1to4Connections_.emplace_back(connect);

        for(int j=0; j<(int)atom1Ptr->bonds_.size(); j++)
        {
            atom2Ptr = atom1Ptr->bonds_[j];
            if(!atom2Ptr || (atom2Ptr == this) || (atom2Ptr == atom1Ptr)) continue;
            connect.atom_ = atom2Ptr;
            connect.numBondsAway_ = 2;
            listOf1to4Connections_.emplace_back(connect);

            for(int k=0; k<(int)atom2Ptr->bonds_.size(); k++)
            {
                atom3Ptr = atom2Ptr->bonds_[k];
                if(!atom3Ptr || (atom3Ptr == this) || (atom3Ptr == atom1Ptr) || (atom3Ptr == atom2Ptr)) continue;
                connect.atom_ = atom3Ptr;
                connect.numBondsAway_ = 3;
                listOf1to4Connections_.emplace_back(connect);

                for(int l=0; l<(int)atom3Ptr->bonds_.size(); l++)
                {
                    atom4Ptr = atom3Ptr->bonds_[l];
                    if(!atom4Ptr || (atom4Ptr == this) || (atom4Ptr == atom1Ptr) || (atom4Ptr == atom2Ptr) || (atom4Ptr == atom3Ptr)) continue;
                    connect.atom_ = atom4Ptr;
                    connect.numBondsAway_ = 4;
                    listOf1to4Connections_.emplace_back(connect);
                }
            }
        }
    }
}

int CAtom::getBondSepTo(const CAtom* atom) const
{
    for(int i=0; i<(int)listOf1to4Connections_.size(); i++)
    {
        if(listOf1to4Connections_[i].atom_ == atom)
            return listOf1to4Connections_[i].numBondsAway_;
    }
    
    return -1;
}

void CAtom::copy(const CAtom& src)
{
    r_ = src.r_;
    Q_ = src.Q_;
    m_ = src.m_;
    sigma_ = src.sigma_;
    ID_ = src.ID_;
    resname_ = src.resname_;
    isMobile_ = src.isMobile_;
    ignoreBondFrom_ = src.ignoreBondFrom_;

    bonds_.clear();
    for(int i=0; i<(int)src.bonds_.size(); i++) bonds_.emplace_back(src.bonds_[i]);

    atomIndex_ = src.atomIndex_;
    molIndex_ = src.molIndex_;
    isSelected_ = src.isSelected_;

    listOf1to4Connections_.clear();
    for(int i=0; i<(int)listOf1to4Connections_.size(); i++) listOf1to4Connections_.emplace_back(src.listOf1to4Connections_[i]);
}

void CAtom::copyIntrinsicAtomProperties(const CAtom &src)
{
    Q_ = src.Q_;
    sigma_ = src.sigma_;
    m_ = src.m_;
    ID_ = src.ID_;
    isMobile_ = src.isMobile_;
    ignoreBondFrom_ = src.ignoreBondFrom_;
}

int CAtom::detectLongestBond(int frame, double& lenghtX, double& lenghtY, double& lenghtZ) const
{
    double dist2, longestDist2 = 0.0;
    int indexLongestBond = -1;
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        dist2 = getDistanceTo2(bonds_[i], frame);
        if(dist2 > longestDist2)
        {
            longestDist2 = dist2;
            indexLongestBond = i;

            lenghtX = fabs(bonds_[i]->r_[frame].x_ - r_[frame].x_);
            lenghtY = fabs(bonds_[i]->r_[frame].y_ - r_[frame].y_);
            lenghtZ = fabs(bonds_[i]->r_[frame].z_ - r_[frame].z_);
        }
    }
    
    return indexLongestBond;
}

void CAtom::findAtomsInMolecule(std::vector<CAtom*>& atomsAtPBCBdry1, std::vector<CAtom*>& atomsAtPBCBdry2, C3DRect pbc, EPBCDir pbcDir, int frame)
{
    int indexCurrBdry = 0;
    double length=0.0;
    double distPBC=0.0, distPBC_2=0.0;
    std::vector<CAtom*>* atomsAtPBCBdry[2] = { &atomsAtPBCBdry1, &atomsAtPBCBdry2 };
    std::map<CAtom*, int> visitedAtoms;
    
    if(pbcDir == dirX) distPBC = pbc.getWidthX();
    if(pbcDir == dirY) distPBC = pbc.getWidthY();
    if(pbcDir == dirZ) distPBC = pbc.getWidthZ();
    distPBC_2 = distPBC / 2.0;
    
    atomsAtPBCBdry1.clear();
    atomsAtPBCBdry2.clear();
    
    visitedAtoms.insert(std::pair<CAtom*, int>(this, 0));
    atomsAtPBCBdry[indexCurrBdry]->emplace_back(this);
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(pbcDir == dirX) length = fabs(bonds_[i]->r_[frame].x_ - r_[frame].x_);
        if(pbcDir == dirY) length = fabs(bonds_[i]->r_[frame].y_ - r_[frame].y_);
        if(pbcDir == dirZ) length = fabs(bonds_[i]->r_[frame].z_ - r_[frame].z_);
        bonds_[i]->findAtomsInMolecule(atomsAtPBCBdry, visitedAtoms, indexCurrBdry, length, pbc, pbcDir, frame);
    }
}

void CAtom::findAtomsInMolecule(std::vector<CAtom*>* atomsAtPBCBdry[], std::map<CAtom*,int>& visitedAtoms, int currBdry, double distFromCaller, C3DRect& pbc, EPBCDir pbcDir, int frame)
{
    double distPBC=0.0, distPBC_2=0.0;
    double length=0.0;
    std::pair<std::map<CAtom*,int>::iterator, bool> ret;

    if(pbcDir == dirX) distPBC = pbc.getWidthX();
    if(pbcDir == dirY) distPBC = pbc.getWidthY();
    if(pbcDir == dirZ) distPBC = pbc.getWidthZ();
    distPBC_2 = distPBC / 2.0;
    
    ret = visitedAtoms.insert(std::pair<CAtom*, int>(this, currBdry));
    if(ret.second == false) return;
    
    if(distFromCaller > distPBC_2) switchBdry(currBdry);
    atomsAtPBCBdry[currBdry]->emplace_back(this);
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        if(pbcDir == dirX) length = fabs(bonds_[i]->r_[frame].x_ - r_[frame].x_);
        if(pbcDir == dirY) length = fabs(bonds_[i]->r_[frame].y_ - r_[frame].y_);
        if(pbcDir == dirZ) length = fabs(bonds_[i]->r_[frame].z_ - r_[frame].z_);
        bonds_[i]->findAtomsInMolecule(atomsAtPBCBdry, visitedAtoms, currBdry, length, pbc, pbcDir, frame);
    }
}

void CAtom::findAtomsInMolecule(std::vector<CAtom*>& atomsInMolecule, int frame)
{
    std::map<CAtom*, int> visitedAtoms;
    std::pair<std::map<CAtom*,int>::iterator, bool> ret;

    atomsInMolecule.clear();
    visitedAtoms.insert(std::pair<CAtom*, int>(this, 0));
    atomsInMolecule.emplace_back(this);

    for(int i=0; i<(int)bonds_.size(); i++)
    {
        bonds_[i]->findAtomsInMolecule(&atomsInMolecule, visitedAtoms, frame);
    }
}

void CAtom::findAtomsInMolecule(std::vector<CAtom*>* atomsInMolecule, std::map<CAtom*,int>& visitedAtoms, int frame)
{
    std::pair<std::map<CAtom*,int>::iterator, bool> ret;

    ret = visitedAtoms.insert(std::pair<CAtom*, int>(this, 0));
    if(ret.second == false) return;

    atomsInMolecule->emplace_back(this);
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        bonds_[i]->findAtomsInMolecule(atomsInMolecule, visitedAtoms, frame);
    }
}

END_CUDA_COMPATIBLE()
