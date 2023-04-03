//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include "MolTwisterStateTools.h"
#include "ProgressBar.h"
#include "MolecularSystemTools.h"
#include "MolecularTools.h"
#include <float.h>

BEGIN_CUDA_COMPATIBLE()

std::shared_ptr<std::vector<int>> CMolTwisterStateTools::getAtomsWithin(double x, double y, double z, double r, C3DRect pbc, const std::vector<CCellListEntry> &atomCellList, int numPBCBisections) const
{
    int cellIndex, numAtoms;
    int iX, iY, iZ;

    auto atomIndices = std::make_shared<std::vector<int>>();
    std::tuple<int, int, int> minIndex = xyzToAtomCellIndex(pbc, x - r, y - r, z - r, numPBCBisections, false);
    std::tuple<int, int, int> maxIndex = xyzToAtomCellIndex(pbc, x + r, y + r, z + r, numPBCBisections, true);

    int minIndexX = std::get<0>(minIndex);
    int minIndexY = std::get<1>(minIndex);
    int minIndexZ = std::get<2>(minIndex);

    int maxIndexX = std::get<0>(maxIndex);
    int maxIndexY = std::get<1>(maxIndex);
    int maxIndexZ = std::get<2>(maxIndex);

    for(int ix=minIndexX; ix<=maxIndexX; ix++)
    {
        for(int iy=minIndexY; iy<=maxIndexY; iy++)
        {
            for(int iz=minIndexZ; iz<=maxIndexZ; iz++)
            {
                // Assign cell inidices
                iX = ix;
                iY = iy;
                iZ = iz;

                // Handle PBCs
                if(ix < 0) iX+= numPBCBisections;
                if(iy < 0) iY+= numPBCBisections;
                if(iz < 0) iZ+= numPBCBisections;

                if(ix >= numPBCBisections) iX-= numPBCBisections;
                if(iy >= numPBCBisections) iY-= numPBCBisections;
                if(iz >= numPBCBisections) iZ-= numPBCBisections;

                if((iX < 0) || (iX >= numPBCBisections)) continue;
                if((iY < 0) || (iY >= numPBCBisections)) continue;
                if((iZ < 0) || (iZ >= numPBCBisections)) continue;

                // Calculate cell index and get atoms from it
                cellIndex = iX + iY*numPBCBisections + iZ*numPBCBisections*numPBCBisections;
                numAtoms = atomCellList[cellIndex].getNumAtoms();
                for(int j=0; j<numAtoms; j++)
                {
                    atomIndices->emplace_back(atomCellList[cellIndex].getAtom(j));
                }
            }
        }
    }

    return atomIndices;
}

std::shared_ptr<std::vector<CMolTwisterStateTools::CCellListEntry>> CMolTwisterStateTools::genAtomCellList(C3DRect pbc, int numPBCBisections, int frame) const
{
    CCellListEntry dummyEntry;
    C3DVector r;
    int size = (numPBCBisections * numPBCBisections * numPBCBisections);
    int cubeIndex;

    auto atomCellList = std::make_shared<std::vector<CCellListEntry>>(size);

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        r = state_->atoms_[i]->r_[frame];
        cubeIndex = xyzToAtomCellIndex(pbc, r.x_, r.y_, r.z_, numPBCBisections);
        if(cubeIndex < 0)
        {
            printf("Waring! Attempt to bisect PBC into smaller blocks failed. Atom found below cube boundaries!\r\n");
            cubeIndex = 0;
        }
        if(cubeIndex >= (int)atomCellList->size())
        {
            printf("Waring! Attempt to bisect PBC into smaller blocks failed. Atom found above cube boundaries!\r\n");
            cubeIndex = (int)atomCellList->size() - 1;
        }

        (*atomCellList)[cubeIndex].addAtom(i);
    }

    return atomCellList;
}

int CMolTwisterStateTools::xyzToAtomCellIndex(C3DRect pbc, double x, double y, double z, int numPBCBisections) const
{
    pbc.expandByFactor(0.1);

    double widthX = pbc.getWidthX();
    double widthY = pbc.getWidthY();
    double widthZ = pbc.getWidthZ();

    int ix = (int)floor(((x - pbc.rLow_.x_) / widthX) * double(numPBCBisections-1));
    int iy = (int)floor(((y - pbc.rLow_.y_) / widthY) * double(numPBCBisections-1));
    int iz = (int)floor(((z - pbc.rLow_.z_) / widthZ) * double(numPBCBisections-1));

    return (ix + iy*numPBCBisections + iz*numPBCBisections*numPBCBisections);
}

std::tuple<int, int, int> CMolTwisterStateTools::xyzToAtomCellIndex(C3DRect pbc, double x, double y, double z, int numPBCBisections, bool getCeiling) const
{
    pbc.expandByFactor(0.1);

    double widthX = pbc.getWidthX();
    double widthY = pbc.getWidthY();
    double widthZ = pbc.getWidthZ();

    int ix, iy, iz;
    if(getCeiling)
    {
        ix = (int)ceil(((x - pbc.rLow_.x_) / widthX) * double(numPBCBisections-1));
        iy = (int)ceil(((y - pbc.rLow_.y_) / widthY) * double(numPBCBisections-1));
        iz = (int)ceil(((z - pbc.rLow_.z_) / widthZ) * double(numPBCBisections-1));
    }
    else
    {
        ix = (int)floor(((x - pbc.rLow_.x_) / widthX) * double(numPBCBisections-1));
        iy = (int)floor(((y - pbc.rLow_.y_) / widthY) * double(numPBCBisections-1));
        iz = (int)floor(((z - pbc.rLow_.z_) / widthZ) * double(numPBCBisections-1));
    }

    return std::tuple<int, int, int>(ix, iy, iz);
}

void CMolTwisterStateTools::generateBonds(double minR, bool verbose, bool, int frame, C3DRect* pbc, std::vector<std::string>* bondAtomsToIgnore)
{
    CProgressBar pb;
    std::vector<int> atomIndices;
    C3DVector R;
    C3DRect PBC;
    double r2 = minR*minR;
    double d2, dWDW1, dWDW2, r2Upper;
    std::string text1;
    std::string text2;
    int numAtomCubeArrayBisections = 50;


    if(!state_ || !state_->view3D_) return;
    if(verbose) printf("\r\n");

    // Note which atoms to ignore bonds from
    if(bondAtomsToIgnore)
    {
        for(int i=0; i<(int)bondAtomsToIgnore->size(); i++)
        {
            for(int j=0; j<(int)state_->atoms_.size(); j++)
            {
                std::string ID = state_->atoms_[j]->getID();
                if(ID == (*bondAtomsToIgnore)[i]) state_->atoms_[j]->ignoreBondFrom_ = true;
            }
        }
    }

    // Delete all existing bonds
    if(frame == -1) frame = state_->currentFrame_;

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        state_->atoms_[i]->deleteAllBonds();
    }

    // We do not wish to loop over all N^2 possible bonds. Therefore, we first create a cube around the
    // PBCs (which is slightly bigger) which we proceed to bisect into several smaller cubes. This cube
    // is stored as an array. We then place each atom within the small boxes (corresponding to physical
    // positions (x,y,z). This is an O(N) algorithm. Thus, we can directly access the group of atoms
    // that are closer than a distance R from the atom we want to perform a bond search on (up to the
    // granularity of the cube bisections). This leasd to an O(N) bond search instead of an O(N^2) search.
    PBC = state_->view3D_->calcPBC(frame);
    std::shared_ptr<std::vector<CMolTwisterStateTools::CCellListEntry>> atomCellList = genAtomCellList(PBC, numAtomCubeArrayBisections, frame);

    // Define largest bond that may occur
    const double largestPossibleBond = 2.32 + 0.05;
    const double largestPossibleBond2 = largestPossibleBond*largestPossibleBond;

    // Create a neighborlist from the atomCellList
    pb.beginProgress("Generate neighborlist");
    std::vector<std::vector<int>> neighborList(state_->atoms_.size());
    for(int i=0,lengthI=(int)state_->atoms_.size(); i<lengthI; i++)
    {
        R = state_->atoms_[i]->r_[frame];
        std::shared_ptr<std::vector<int>> atomIndices = getAtomsWithin(R.x_, R.y_, R.z_, 6.0, PBC, *atomCellList, numAtomCubeArrayBisections);
        neighborList[i].reserve(atomIndices->size());
        for(size_t j=0,lengthJ=(int)atomIndices->size(); j<lengthJ; j++)
        {
            int atomIndex = (*atomIndices)[j];
            C3DVector r = state_->atoms_[atomIndex]->r_[frame];
            if(R.distToAcrossPBC2(r, PBC) <= largestPossibleBond2)
            {
                neighborList[i].emplace_back(atomIndex);
            }
        }
        pb.updateProgress(i, (int)state_->atoms_.size());
    }
    pb.endProgress();

    // Part 1 of [Journal of Cheminformatics 2012, 4,26] algorithm : generate
    // a bond if dMinR < r_ij < r_i+r_j+0.4, where r_i and r_j are covalent
    // radii and r_ij is the distance between atom i and atom j
    pb.beginProgress("Bond gen step 1/3");
    for(int i=0,lengthI=(int)state_->atoms_.size(); i<lengthI; i++)
    {
        if(state_->atoms_[i]->ignoreBondFrom_) continue;

        for(int k=0,lengthK=(int)neighborList[i].size(); k<lengthK; k++)
        {
            int j = neighborList[i][k];
            if(i == j) continue;

            if(state_->atoms_[j]->ignoreBondFrom_) continue;

            if(!pbc) d2 = state_->atoms_[i]->getDistanceTo2(state_->atoms_[j].get(), frame);
            else     d2 = state_->atoms_[i]->getDistanceTo2UsingPBC(state_->atoms_[j].get(), frame, *pbc);

            if(d2 > r2)
            {
                text1 = state_->atoms_[i]->getID();
                text2 = state_->atoms_[j]->getID();

                dWDW1 = state_->defaultAtProp_.getCovalentRadius(text1.data());
                dWDW2 = state_->defaultAtProp_.getCovalentRadius(text2.data());
                r2Upper = dWDW1 + dWDW2 + 0.4;
                r2Upper*= r2Upper;

                if(d2 < r2Upper)
                {
                    state_->atoms_[i]->attachBondTo(state_->atoms_[j].get());
                    if(verbose)
                        printf("\tBond detected: %s[%i] - %s[%i]\r\n", text1.data(), i, text2.data(), j);
                }
            }
        }
        pb.updateProgress(i, (int)state_->atoms_.size());
    }
    pb.endProgress();


    // Part 2 of [Journal of Cheminformatics 2012, 4,26] algorithm : if number of bonds
    // connected to C, N, P or S have more than 4 bonds connected, the longest bonds
    // are removed until only 4 bonds are left
    double lengthX, lengthY, lengthZ;
    int locatedAtomIndex;
    int indexOfLongestBond;
    int atomIndexC = state_->defaultAtProp_.identifyAtom("C");
    int atomIndexN = state_->defaultAtProp_.identifyAtom("N");
    int atomIndexP = state_->defaultAtProp_.identifyAtom("P");
    int atomIndexS = state_->defaultAtProp_.identifyAtom("S");
    pb.beginProgress("Bond gen step 2/3");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        std::string  ID = state_->atoms_[i]->getID();
        locatedAtomIndex = state_->defaultAtProp_.identifyAtom(ID.data());

        if((locatedAtomIndex == atomIndexC) || (locatedAtomIndex == atomIndexN) ||
           (locatedAtomIndex == atomIndexP) || (locatedAtomIndex == atomIndexS))
        {
            while(state_->atoms_[i]->getNumBonds() > 4)
            {
                indexOfLongestBond = state_->atoms_[i]->detectLongestBond(frame, lengthX, lengthY, lengthZ);
                CAtom* atom = state_->atoms_[i]->getBondDest(indexOfLongestBond);
                if(atom)
                {
                    state_->atoms_[i]->deleteBondTo(atom);
                    atom->deleteBondTo(state_->atoms_[i].get());
                }
            }
        }
        pb.updateProgress(i, (int)state_->atoms_.size());
    }
    pb.endProgress();

    state_->genMolIndices();

    pb.beginProgress("Bond gen step 3/3");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]) state_->atoms_[i]->buildListOf1to4Connections();
        pb.updateProgress(i, (int)state_->atoms_.size());
    }
    pb.endProgress();
}

void CMolTwisterStateTools::getAtomsOfMolecule(int molIndex, std::vector<int>& atomIndices) const
{
    atomIndices.clear();
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]->getMolIndex() == molIndex)
            atomIndices.emplace_back(i);
    }
}

C3DVector CMolTwisterStateTools::getMoleculeDipoleMoment(const std::vector<int>& atomIndices, int frameIndex, C3DVector& Rc, bool chargeNeutralFormulation) const
{
    C3DVector P, Rc_int;

    Rc = getCenterOfMass(atomIndices, frameIndex);
    if(!chargeNeutralFormulation) Rc_int = Rc;

    for(int i=0; i<(int)atomIndices.size(); i++)
    {
        double q_i = state_->atoms_[atomIndices[i]]->Q_;
        C3DVector r_i = state_->atoms_[atomIndices[i]]->r_[frameIndex];
        C3DVector d_i = r_i - Rc_int;
        P+= (d_i * q_i);
    }

    return P;
}

int CMolTwisterStateTools::getNumMolecules() const
{
    int largestIndex = -1;

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]->getMolIndex() > largestIndex)
            largestIndex = state_->atoms_[i]->getMolIndex();
    }

    return (largestIndex + 1);
}

void CMolTwisterStateTools::getMolecules(std::vector<std::vector<int>>& molecules, const std::vector<std::string>& atomTypesToInclude) const
{
    std::map<std::string, int> mapIDtoIndex;

    // Prepare a map that can be used to check if an atom type is to be included or not
    int numAtomTypesToInclude = (int)atomTypesToInclude.size();
    for(int i=0; i<numAtomTypesToInclude; i++) mapIDtoIndex[atomTypesToInclude[i]] = 1;

    // Loop through all molecules and build a list, excluding all atom types not present
    // in aAtomTypesToInclude, or include all is size(aAtomTypesToInclude) == 0.
    int numMolecules = getNumMolecules();
    for(int i=0; i<numMolecules; i++)
    {
        // Get atom indices in molecule i
        std::vector<int> atomIndices;
        getAtomsOfMolecule(i, atomIndices);
        if(atomIndices.size() <= 0) continue;

        // Exclude the indices not belonging to aAtomTypesToInclude (or include all)
        if(numAtomTypesToInclude != 0)
        {
            for(int j=0; j<(int)atomIndices.size(); j++)
            {
                std::string ID = state_->atoms_[atomIndices[j]]->getID();
                if(mapIDtoIndex.find(ID) == mapIDtoIndex.end())
                {
                    atomIndices.erase(atomIndices.begin() + j);
                    j--;
                }
            }
        }

        // Include molecule into list, if there are any indices left in the molecule
        if(atomIndices.size() > 0)
        {
            molecules.emplace_back(atomIndices);
        }
    }
}

void CMolTwisterStateTools::atomicUnwrap(C3DRect pbc)
{
    int numAtomIndices, numFrames, numMolecules = -1;
    double PBCx = pbc.getWidthX();
    double PBCy = pbc.getWidthY();
    double PBCz = pbc.getWidthZ();
    std::vector<int> atomIndices;
    std::vector<CAtom*> atomsAtPBCBdry1, atomsAtPBCBdry2;


    if(state_->atoms_.size() < 1) return;
    numFrames = (int)state_->atoms_[0]->r_.size();


    // Count number of molecules
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]->getMolIndex() > numMolecules)
            numMolecules = state_->atoms_[i]->getMolIndex();
    }

    numMolecules++;


    // Unwrap the system
    for(int iFrame=0; iFrame<numFrames; iFrame++)
    {
        for(int iMol=0; iMol<numMolecules; iMol++)
        {
            getAtomsOfMolecule(iMol, atomIndices);
            numAtomIndices = (int)atomIndices.size();
            if(numAtomIndices <= 0) continue;


            // Unwrap x-coordinates of molecule : 1) find one set of atoms close to one side of the boundary and another set that
            // is close to the other side of the boundary 2) find how far each of the sets span along x inside the PBCs [i.e. max
            // X, min X] 3) if the high value of span A is lower than the low value of span B, then set A of atoms should be shifted
            // to the higher x-boundary, else set B should be moved to the higher x-boundary
            state_->atoms_[atomIndices[0]]->findAtomsInMolecule(atomsAtPBCBdry1, atomsAtPBCBdry2, pbc, CAtom::dirX, iFrame);

            double boundsOfAtomsAtPBCBdry1[2] = { DBL_MAX, -DBL_MAX };
            double boundsOfAtomsAtPBCBdry2[2] = { DBL_MAX, -DBL_MAX };
            for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++)
            {
                if(atomsAtPBCBdry1[i]->r_[iFrame].x_ < boundsOfAtomsAtPBCBdry1[0]) boundsOfAtomsAtPBCBdry1[0] = atomsAtPBCBdry1[i]->r_[iFrame].x_;
                if(atomsAtPBCBdry1[i]->r_[iFrame].x_ > boundsOfAtomsAtPBCBdry1[1]) boundsOfAtomsAtPBCBdry1[1] = atomsAtPBCBdry1[i]->r_[iFrame].x_;
            }
            for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++)
            {
                if(atomsAtPBCBdry2[i]->r_[iFrame].x_ < boundsOfAtomsAtPBCBdry2[0]) boundsOfAtomsAtPBCBdry2[0] = atomsAtPBCBdry2[i]->r_[iFrame].x_;
                if(atomsAtPBCBdry2[i]->r_[iFrame].x_ > boundsOfAtomsAtPBCBdry2[1]) boundsOfAtomsAtPBCBdry2[1] = atomsAtPBCBdry2[i]->r_[iFrame].x_;
            }

            if((atomsAtPBCBdry1.size() != 0) && (atomsAtPBCBdry2.size() != 0))
            {
                if(boundsOfAtomsAtPBCBdry1[1] < boundsOfAtomsAtPBCBdry2[0])
                {
                    for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++) atomsAtPBCBdry1[i]->r_[iFrame].x_+= PBCx;
                }
                else
                {
                    for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++) atomsAtPBCBdry2[i]->r_[iFrame].x_+= PBCx;
                }
            }


            // Unwrap y-coordinates of molecule
            state_->atoms_[atomIndices[0]]->findAtomsInMolecule(atomsAtPBCBdry1, atomsAtPBCBdry2, pbc, CAtom::dirY, iFrame);

            boundsOfAtomsAtPBCBdry1[0] = DBL_MAX; boundsOfAtomsAtPBCBdry1[1] = -DBL_MAX;
            boundsOfAtomsAtPBCBdry2[0] = DBL_MAX; boundsOfAtomsAtPBCBdry2[1] = -DBL_MAX;
            for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++)
            {
                if(atomsAtPBCBdry1[i]->r_[iFrame].y_ < boundsOfAtomsAtPBCBdry1[0]) boundsOfAtomsAtPBCBdry1[0] = atomsAtPBCBdry1[i]->r_[iFrame].y_;
                if(atomsAtPBCBdry1[i]->r_[iFrame].y_ > boundsOfAtomsAtPBCBdry1[1]) boundsOfAtomsAtPBCBdry1[1] = atomsAtPBCBdry1[i]->r_[iFrame].y_;
            }
            for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++)
            {
                if(atomsAtPBCBdry2[i]->r_[iFrame].y_ < boundsOfAtomsAtPBCBdry2[0]) boundsOfAtomsAtPBCBdry2[0] = atomsAtPBCBdry2[i]->r_[iFrame].y_;
                if(atomsAtPBCBdry2[i]->r_[iFrame].y_ > boundsOfAtomsAtPBCBdry2[1]) boundsOfAtomsAtPBCBdry2[1] = atomsAtPBCBdry2[i]->r_[iFrame].y_;
            }

            if(boundsOfAtomsAtPBCBdry1[1] < boundsOfAtomsAtPBCBdry2[0])
            {
                for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++) atomsAtPBCBdry1[i]->r_[iFrame].y_+= PBCy;
            }
            else
            {
                for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++) atomsAtPBCBdry2[i]->r_[iFrame].y_+= PBCy;
            }


            // Unwrap z-coordinates of molecule
            state_->atoms_[atomIndices[0]]->findAtomsInMolecule(atomsAtPBCBdry1, atomsAtPBCBdry2, pbc, CAtom::dirZ, iFrame);

            boundsOfAtomsAtPBCBdry1[0] = DBL_MAX; boundsOfAtomsAtPBCBdry1[1] = -DBL_MAX;
            boundsOfAtomsAtPBCBdry2[0] = DBL_MAX; boundsOfAtomsAtPBCBdry2[1] = -DBL_MAX;
            for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++)
            {
                if(atomsAtPBCBdry1[i]->r_[iFrame].z_ < boundsOfAtomsAtPBCBdry1[0]) boundsOfAtomsAtPBCBdry1[0] = atomsAtPBCBdry1[i]->r_[iFrame].z_;
                if(atomsAtPBCBdry1[i]->r_[iFrame].z_ > boundsOfAtomsAtPBCBdry1[1]) boundsOfAtomsAtPBCBdry1[1] = atomsAtPBCBdry1[i]->r_[iFrame].z_;
            }
            for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++)
            {
                if(atomsAtPBCBdry2[i]->r_[iFrame].z_ < boundsOfAtomsAtPBCBdry2[0]) boundsOfAtomsAtPBCBdry2[0] = atomsAtPBCBdry2[i]->r_[iFrame].z_;
                if(atomsAtPBCBdry2[i]->r_[iFrame].z_ > boundsOfAtomsAtPBCBdry2[1]) boundsOfAtomsAtPBCBdry2[1] = atomsAtPBCBdry2[i]->r_[iFrame].z_;
            }

            if(boundsOfAtomsAtPBCBdry1[1] < boundsOfAtomsAtPBCBdry2[0])
            {
                for(int i=0; i<(int)atomsAtPBCBdry1.size(); i++) atomsAtPBCBdry1[i]->r_[iFrame].z_+= PBCz;
            }
            else
            {
                for(int i=0; i<(int)atomsAtPBCBdry2.size(); i++) atomsAtPBCBdry2[i]->r_[iFrame].z_+= PBCz;
            }
        }
    }
}

void CMolTwisterStateTools::getAllVisibleBondsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2) const
{
    std::map<CAtom*,int> mapAtomPtrToIndex;
    CAtom* atom1;
    CAtom* atom2;
    int numBonds1;
    int indexAtom1, indexAtom2;

    atoms1.clear();
    atoms2.clear();

    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atom1 = state_->atoms_[i].get();
        numBonds1 = atom1->getNumBonds();
        for(int j=0; j<numBonds1; j++)
        {
            atom2 = atom1->getBondDest(j);

            indexAtom1 = mapAtomPtrToIndex[atom1];
            indexAtom2 = mapAtomPtrToIndex[atom2];

            atoms1.emplace_back(indexAtom1);
            atoms2.emplace_back(indexAtom2);
        }
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateBondIndices(atoms1, atoms2);
}

void CMolTwisterStateTools::getAllMDBondsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC)
{
    CMolecularTools molTools(state_, stdOut_);
    CProgressBar pb;
    std::map<CAtom*,int> mapAtomPtrToIndex;
    C3DRect PBC;
    std::string stringAtom1, stringAtom2;
    double R0, R02;
    CAtom* atom1;
    CAtom* atom2;
    char progbarHeader[128];
    int indexAtom1, indexAtom2;

    atoms1.clear();
    atoms2.clear();
    mdTypeIndex.clear();

    if(state_->currentFrame_ < 0) return;
    if(!state_->view3D_) return;
    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);
    PBC = state_->view3D_->calcPBC(state_->currentFrame_);

    int numBondSpecs = state_->mdFFBondList_.size();
    for(int typeIndex=0; typeIndex<numBondSpecs; typeIndex++)
    {
        CMDFFBond* bondType = state_->mdFFBondList_.get(typeIndex);
        if(!bondType) continue;

        stringAtom1 = bondType->getAtomInBond(0);
        stringAtom2 = bondType->getAtomInBond(1);

        sprintf(progbarHeader, "Bond search %i/%i (%s-%s)", typeIndex+1, numBondSpecs, stringAtom1.data(), stringAtom2.data());
        pb.beginProgress(progbarHeader);

        CMDFFBond::EDetCrit detCrit = bondType->getBondDetectionCriteria(R0); R02 = R0*R0;
        if(detCrit == CMDFFBond::critAllRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsType2;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            state_->getAtomsWithID(stringAtom2, atomsType2);

            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                for(int k=0; k<(int)atomsType2.size(); k++)
                {
                    if(atomsType1[j] == atomsType2[k]) continue;
                    if(molTools.calcDistance2(atomsType1[j], atomsType2[k], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                    {
                        indexAtom1 = mapAtomPtrToIndex[atomsType1[j]];
                        indexAtom2 = mapAtomPtrToIndex[atomsType2[k]];

                        atoms1.emplace_back(indexAtom1);
                        atoms2.emplace_back(indexAtom2);
                        mdTypeIndex.emplace_back(typeIndex);
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFBond::critMolRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsInMolecule1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                atom1->findAtomsInMolecule(atomsInMolecule1, state_->currentFrame_);
                for(int k=0; k<(int)atomsInMolecule1.size(); k++)
                {
                    atom2 = atomsInMolecule1[k];
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(atom1 == atom2) continue;
                    if((ID == stringAtom2) && (molTools.calcDistance2(atom1, atom2, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                    {
                        indexAtom1 = mapAtomPtrToIndex[atom1];
                        indexAtom2 = mapAtomPtrToIndex[atom2];

                        atoms1.emplace_back(indexAtom1);
                        atoms2.emplace_back(indexAtom2);
                        mdTypeIndex.emplace_back(typeIndex);
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFBond::critOnlyVisibleBonds)
        {
            std::vector<CAtom*> atomsType1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                for(int k=0; k<atom1->getNumBonds(); k++)
                {
                    atom2 = atom1->getBondDest(k);
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(ID == stringAtom2)
                    {
                        indexAtom1 = mapAtomPtrToIndex[atom1];
                        indexAtom2 = mapAtomPtrToIndex[atom2];

                        atoms1.emplace_back(indexAtom1);
                        atoms2.emplace_back(indexAtom2);
                        mdTypeIndex.emplace_back(typeIndex);
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFBond::critOnly14Bonds)
        {
            std::vector<CAtom*> atomsType1;
            CAtom* atom3;
            CAtom* atom4;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                for(int k=0; k<atom1->getNumBonds(); k++)
                {
                    atom2 = atom1->getBondDest(k);
                    if(!atom2) continue;

                    for(int l=0; l<atom2->getNumBonds(); l++)
                    {
                        atom3 = atom2->getBondDest(l);
                        if(!atom3) continue;

                        for(int m=0; m<atom3->getNumBonds(); m++)
                        {
                            atom4 = atom3->getBondDest(m);
                            if(!atom4) continue;

                            std::string ID = atom4->getID();
                            if(ID == stringAtom2)
                            {
                                indexAtom1 = mapAtomPtrToIndex[atom1];
                                indexAtom2 = mapAtomPtrToIndex[atom4];

                                if((atom4 != atom3) && (atom4 != atom2) && (atom4 != atom1))
                                {
                                    atoms1.emplace_back(indexAtom1);
                                    atoms2.emplace_back(indexAtom2);
                                    mdTypeIndex.emplace_back(typeIndex);
                                }
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else
        {
            printf("Error occured in bond search: unsupported determination criteria used in bond search!\r\n");
            return;
        }

        pb.endProgress();
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateBondIndices(atoms1, atoms2, &mdTypeIndex);
}

void CMolTwisterStateTools::getAllVisibleAnglesInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3) const
{
    std::map<CAtom*,int> mapAtomPtrToIndex;
    CAtom* atom1;
    CAtom* atom2;
    CAtom* atom3;
    int numBonds1, numBonds2;
    int indexAtom1, indexAtom2, indexAtom3;

    atoms1.clear();
    atoms2.clear();
    atoms3.clear();

    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atom1 = state_->atoms_[i].get();
        numBonds1 = atom1->getNumBonds();
        for(int j=0; j<numBonds1; j++)
        {
            atom2 = atom1->getBondDest(j);
            numBonds2 = atom2->getNumBonds();
            for(int k=0; k<numBonds2; k++)
            {
                atom3 = atom2->getBondDest(k);

                indexAtom1 = mapAtomPtrToIndex[atom1];
                indexAtom2 = mapAtomPtrToIndex[atom2];
                indexAtom3 = mapAtomPtrToIndex[atom3];

                atoms1.emplace_back(indexAtom1);
                atoms2.emplace_back(indexAtom2);
                atoms3.emplace_back(indexAtom3);
            }
        }
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateAngleIndices(atoms1, atoms2, atoms3);
}

void CMolTwisterStateTools::getAllMDAnglesInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC)
{
    CMolecularTools molTools(state_, stdOut_);
    CProgressBar pb;
    std::map<CAtom*,int> mapAtomPtrToIndex;
    C3DRect PBC;
    std::string stringAtom1, stringAtom2, stringAtom3;
    CAtom* atom1;
    CAtom* atom2;
    CAtom* atom3;
    double R0, R02;
    char progbarHeader[128];
    int indexAtom1, indexAtom2, indexAtom3;

    atoms1.clear();
    atoms2.clear();
    atoms3.clear();
    mdTypeIndex.clear();

    if(state_->currentFrame_ < 0) return;
    if(!state_->view3D_) return;
    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);
    PBC = state_->view3D_->calcPBC(state_->currentFrame_);

    int numAngleSpecs = state_->mdFFAngleList_.size();
    for(int typeIndex=0; typeIndex<numAngleSpecs; typeIndex++)
    {
        CMDFFAngle* angleType = state_->mdFFAngleList_.get(typeIndex);
        if(!angleType) continue;

        stringAtom1 = angleType->getAtomInBond(0);
        stringAtom2 = angleType->getAtomInBond(1);
        stringAtom3 = angleType->getAtomInBond(2);

        sprintf(progbarHeader, "Angle search %i/%i (%s-%s-%s)", typeIndex+1, numAngleSpecs, stringAtom1.data(), stringAtom2.data(), stringAtom3.data());
        pb.beginProgress(progbarHeader);

        CMDFFAngle::EDetCrit detCrit = angleType->getBondDetectionCriteria(R0); R02 = R0*R0;

        if(detCrit == CMDFFAngle::critAllRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsType2, atomsType3;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            state_->getAtomsWithID(stringAtom2, atomsType2);
            state_->getAtomsWithID(stringAtom3, atomsType3);

            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                for(int k=0; k<(int)atomsType2.size(); k++)
                {
                    if(atomsType1[j] == atomsType2[k]) continue;
                    if(molTools.calcDistance2(atomsType1[j], atomsType2[k], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                    {
                        for(int l=0; l<(int)atomsType3.size(); l++)
                        {
                            if(atomsType3[l] == atomsType2[k]) continue;
                            if(atomsType3[l] == atomsType1[j]) continue;
                            if(molTools.calcDistance2(atomsType2[k], atomsType3[l], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                            {
                                indexAtom1 = mapAtomPtrToIndex[atomsType1[j]];
                                indexAtom2 = mapAtomPtrToIndex[atomsType2[k]];
                                indexAtom3 = mapAtomPtrToIndex[atomsType3[l]];

                                atoms1.emplace_back(indexAtom1);
                                atoms2.emplace_back(indexAtom2);
                                atoms3.emplace_back(indexAtom3);
                                mdTypeIndex.emplace_back(typeIndex);
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFAngle::critMolRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsInMolecule1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                atom1->findAtomsInMolecule(atomsInMolecule1, state_->currentFrame_);
                for(int k=0; k<(int)atomsInMolecule1.size(); k++)
                {
                    atom2 = atomsInMolecule1[k];
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(atom1 == atom2) continue;
                    if((ID == stringAtom2) && (molTools.calcDistance2(atom1, atom2, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                    {
                        for(int l=0; l<(int)atomsInMolecule1.size(); l++)
                        {
                            atom3 = atomsInMolecule1[l];
                            if(!atom3) continue;

                            ID = atom3->getID();
                            if(atom3 == atom1) continue;
                            if(atom3 == atom2) continue;
                            if((ID == stringAtom3) && (molTools.calcDistance2(atom2, atom3, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                            {
                                indexAtom1 = mapAtomPtrToIndex[atom1];
                                indexAtom2 = mapAtomPtrToIndex[atom2];
                                indexAtom3 = mapAtomPtrToIndex[atom3];

                                atoms1.emplace_back(indexAtom1);
                                atoms2.emplace_back(indexAtom2);
                                atoms3.emplace_back(indexAtom3);
                                mdTypeIndex.emplace_back(typeIndex);
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFAngle::critOnlyVisibleBonds)
        {
            std::vector<CAtom*> atomsType1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                for(int k=0; k<atom1->getNumBonds(); k++)
                {
                    atom2 = atom1->getBondDest(k);
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(ID == stringAtom2)
                    {
                        for(int l=0; l<atom2->getNumBonds(); l++)
                        {
                            atom3 = atom2->getBondDest(l);
                            if(!atom3) continue;

                            ID = atom3->getID();
                            if(ID == stringAtom3)
                            {
                                indexAtom1 = mapAtomPtrToIndex[atom1];
                                indexAtom2 = mapAtomPtrToIndex[atom2];
                                indexAtom3 = mapAtomPtrToIndex[atom3];

                                atoms1.emplace_back(indexAtom1);
                                atoms2.emplace_back(indexAtom2);
                                atoms3.emplace_back(indexAtom3);
                                mdTypeIndex.emplace_back(typeIndex);
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else
        {
            printf("Error occured in angle search: unsupported determination criteria used in bond search!\r\n");
            return;
        }

        pb.endProgress();
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateAngleIndices(atoms1, atoms2, atoms3, &mdTypeIndex);
}

void CMolTwisterStateTools::getAllVisibleDihedralsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4) const
{
    std::map<CAtom*,int> mapAtomPtrToIndex;
    CAtom* atom1;
    CAtom* atom2;
    CAtom* atom3;
    CAtom* atom4;
    int numBonds1, numBonds2, numBonds3;
    int indexAtom1, indexAtom2, indexAtom3, indexAtom4;

    atoms1.clear();
    atoms2.clear();
    atoms3.clear();
    atoms4.clear();

    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);

    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atom1 = state_->atoms_[i].get();
        numBonds1 = atom1->getNumBonds();
        for(int j=0; j<numBonds1; j++)
        {
            atom2 = atom1->getBondDest(j);
            numBonds2 = atom2->getNumBonds();
            for(int k=0; k<numBonds2; k++)
            {
                atom3 = atom2->getBondDest(k);
                numBonds3 = atom3->getNumBonds();
                for(int l=0; l<numBonds3; l++)
                {
                    atom4 = atom3->getBondDest(l);

                    indexAtom1 = mapAtomPtrToIndex[atom1];
                    indexAtom2 = mapAtomPtrToIndex[atom2];
                    indexAtom3 = mapAtomPtrToIndex[atom3];
                    indexAtom4 = mapAtomPtrToIndex[atom4];

                    atoms1.emplace_back(indexAtom1);
                    atoms2.emplace_back(indexAtom2);
                    atoms3.emplace_back(indexAtom3);
                    atoms4.emplace_back(indexAtom4);
                }
            }
        }
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateDihIndices(atoms1, atoms2, atoms3, atoms4);
}

void CMolTwisterStateTools::getAllMDDihedralsInSystem(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4, std::vector<int>& mdTypeIndex, bool bondsAcrossPBC)
{
    CMolecularTools molTools(state_, stdOut_);
    CProgressBar pb;
    std::map<CAtom*,int> mapAtomPtrToIndex;
    C3DRect PBC;
    std::string stringAtom1, stringAtom2, stringAtom3, stringAtom4;
    CAtom* atom1;
    CAtom* atom2;
    CAtom* atom3;
    CAtom* atom4;
    double R0, R02;
    char progbarHeader[128];
    int indexAtom1, indexAtom2, indexAtom3, indexAtom4;

    atoms1.clear();
    atoms2.clear();
    atoms3.clear();
    atoms4.clear();
    mdTypeIndex.clear();

    if(state_->currentFrame_ < 0) return;
    if(!state_->view3D_) return;
    state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);
    PBC = state_->view3D_->calcPBC(state_->currentFrame_);

    int iNumDihSpecs = state_->mdFFDihList_.size();
    for(int iTypeIndex=0; iTypeIndex<iNumDihSpecs; iTypeIndex++)
    {
        CMDFFDih* dihType = state_->mdFFDihList_.get(iTypeIndex);
        if(!dihType) continue;

        stringAtom1 = dihType->getAtomInBond(0);
        stringAtom2 = dihType->getAtomInBond(1);
        stringAtom3 = dihType->getAtomInBond(2);
        stringAtom4 = dihType->getAtomInBond(3);

        sprintf(progbarHeader, "Dihedral search %i/%i (%s-%s-%s-%s)", iTypeIndex+1, iNumDihSpecs, stringAtom1.data(), stringAtom2.data(), stringAtom3.data(), stringAtom4.data());
        pb.beginProgress(progbarHeader);

        CMDFFDih::EDetCrit detCrit = dihType->getBondDetectionCriteria(R0); R02 = R0*R0;

        if(detCrit == CMDFFDih::critAllRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsType2, atomsType3, atomsType4;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            state_->getAtomsWithID(stringAtom2, atomsType2);
            state_->getAtomsWithID(stringAtom3, atomsType3);
            state_->getAtomsWithID(stringAtom4, atomsType4);

            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                for(int k=0; k<(int)atomsType2.size(); k++)
                {
                    if(atomsType1[j] == atomsType2[k]) continue;
                    if(molTools.calcDistance2(atomsType1[j], atomsType2[k], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                    {
                        for(int l=0; l<(int)atomsType3.size(); l++)
                        {
                            if(atomsType3[l] == atomsType2[k]) continue;
                            if(atomsType3[l] == atomsType1[j]) continue;
                            if(molTools.calcDistance2(atomsType2[k], atomsType3[l], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                            {
                                for(int m=0; m<(int)atomsType4.size(); m++)
                                {
                                    if(atomsType4[m] == atomsType3[l]) continue;
                                    if(atomsType4[m] == atomsType2[k]) continue;
                                    if(atomsType4[m] == atomsType1[j]) continue;
                                    if(molTools.calcDistance2(atomsType3[l], atomsType4[m], state_->currentFrame_, PBC, bondsAcrossPBC) < R02)
                                    {
                                        indexAtom1 = mapAtomPtrToIndex[atomsType1[j]];
                                        indexAtom2 = mapAtomPtrToIndex[atomsType2[k]];
                                        indexAtom3 = mapAtomPtrToIndex[atomsType3[l]];
                                        indexAtom4 = mapAtomPtrToIndex[atomsType4[m]];

                                        atoms1.emplace_back(indexAtom1);
                                        atoms2.emplace_back(indexAtom2);
                                        atoms3.emplace_back(indexAtom3);
                                        atoms4.emplace_back(indexAtom4);
                                        mdTypeIndex.emplace_back(iTypeIndex);
                                    }
                                }
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFDih::critMolRLessThan_R0)
        {
            std::vector<CAtom*> atomsType1, atomsInMolecule1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                atom1->findAtomsInMolecule(atomsInMolecule1, state_->currentFrame_);
                for(int k=0; k<(int)atomsInMolecule1.size(); k++)
                {
                    atom2 = atomsInMolecule1[k];
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(atom1 == atom2) continue;
                    if((ID == stringAtom2) && (molTools.calcDistance2(atom1, atom2, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                    {
                        for(int l=0; l<(int)atomsInMolecule1.size(); l++)
                        {
                            atom3 = atomsInMolecule1[l];
                            if(!atom3) continue;

                            ID = atom3->getID();
                            if(atom3 == atom1) continue;
                            if(atom3 == atom2) continue;
                            if((ID == stringAtom3) && (molTools.calcDistance2(atom2, atom3, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                            {
                                for(int m=0; m<(int)atomsInMolecule1.size(); m++)
                                {
                                    atom4 = atomsInMolecule1[m];
                                    if(!atom4) continue;

                                    ID = atom4->getID();
                                    if(atom4 == atom1) continue;
                                    if(atom4 == atom2) continue;
                                    if(atom4 == atom3) continue;
                                    if((ID == stringAtom4) && (molTools.calcDistance2(atom3, atom4, state_->currentFrame_, PBC, bondsAcrossPBC) < R02))
                                    {
                                        indexAtom1 = mapAtomPtrToIndex[atom1];
                                        indexAtom2 = mapAtomPtrToIndex[atom2];
                                        indexAtom3 = mapAtomPtrToIndex[atom3];
                                        indexAtom4 = mapAtomPtrToIndex[atom4];

                                        atoms1.emplace_back(indexAtom1);
                                        atoms2.emplace_back(indexAtom2);
                                        atoms3.emplace_back(indexAtom3);
                                        atoms4.emplace_back(indexAtom4);
                                        mdTypeIndex.emplace_back(iTypeIndex);
                                    }
                                }
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else if(detCrit == CMDFFDih::critOnlyVisibleBonds)
        {
            std::vector<CAtom*> atomsType1;

            state_->getAtomsWithID(stringAtom1, atomsType1);
            for(int j=0; j<(int)atomsType1.size(); j++)
            {
                atom1 = atomsType1[j];
                if(!atom1) continue;

                for(int k=0; k<atom1->getNumBonds(); k++)
                {
                    atom2 = atom1->getBondDest(k);
                    if(!atom2) continue;

                    std::string ID = atom2->getID();
                    if(ID == stringAtom2)
                    {
                        for(int l=0; l<atom2->getNumBonds(); l++)
                        {
                            atom3 = atom2->getBondDest(l);
                            if(!atom3) continue;

                            ID = atom3->getID();
                            if(ID == stringAtom3)
                            {
                                for(int m=0; m<atom3->getNumBonds(); m++)
                                {
                                    atom4 = atom3->getBondDest(m);
                                    if(!atom4) continue;

                                    ID = atom4->getID();
                                    if(ID == stringAtom4)
                                    {
                                        indexAtom1 = mapAtomPtrToIndex[atom1];
                                        indexAtom2 = mapAtomPtrToIndex[atom2];
                                        indexAtom3 = mapAtomPtrToIndex[atom3];
                                        indexAtom4 = mapAtomPtrToIndex[atom4];

                                        atoms1.emplace_back(indexAtom1);
                                        atoms2.emplace_back(indexAtom2);
                                        atoms3.emplace_back(indexAtom3);
                                        atoms4.emplace_back(indexAtom4);
                                        mdTypeIndex.emplace_back(iTypeIndex);
                                    }
                                }
                            }
                        }
                    }
                }

                pb.updateProgress(j, (int)atomsType1.size());
            }
        }
        else
        {
            printf("Error occured in dihedral search: unsupported determination criteria used in bond search!\r\n");
            return;
        }

        pb.endProgress();
    }

    CMolecularSystemTools(state_, stdOut_).removeDuplicateDihIndices(atoms1, atoms2, atoms3, atoms4, &mdTypeIndex);
}

void CMolTwisterStateTools::reportConsistencyOfMDForceField(FILE* file)
{
    std::vector<int> bondAtoms1, bondAtoms2, bondMDTypeIndices;
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    std::vector<int> dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices;


    // Find all bonds, angles and dihedrals
    getAllMDBondsInSystem(bondAtoms1, bondAtoms2, bondMDTypeIndices);
    getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices);
    getAllMDDihedralsInSystem(dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices);


    // Check inconsistencies between bonds found in bond-search and
    // link to bond types. All bond types should be in use and no
    // non-existent bond types should be used.
    std::vector<bool> bondTypesUsedByDataFile;
    bondTypesUsedByDataFile.resize(state_->mdFFBondList_.size(), false);

    for(int i=0; i<(int)bondAtoms1.size(); i++)
    {
        if(bondMDTypeIndices[i] < (int)bondTypesUsedByDataFile.size()) bondTypesUsedByDataFile[bondMDTypeIndices[i]] = true;
        else
        {
            fprintf(file, "\tBond type index %i does not exist, there are only %i bond types declared in total!\r\n", bondMDTypeIndices[i]+1, (int)bondTypesUsedByDataFile.size());
        }
    }


    // Check inconsistencies between angles found in angle-search and
    // link to angle types. All angle types should be in use and no
    // non-existent angle types should be used.
    std::vector<bool> angleTypesUsedByDataFile;
    angleTypesUsedByDataFile.resize(state_->mdFFAngleList_.size(), false);

    fprintf(file, "\r\n");

    for(int i=0; i<(int)angleAtoms1.size(); i++)
    {
        if(angleMDTypeIndices[i] < (int)angleTypesUsedByDataFile.size()) angleTypesUsedByDataFile[angleMDTypeIndices[i]] = true;
        else
        {
            fprintf(file, "\tAngle type index %i does not exist, there are only %i angle types declared in total!\r\n", angleMDTypeIndices[i]+1, (int)angleTypesUsedByDataFile.size());
        }
    }


    // Check inconsistencies between dihedrals found in dihedral-search and
    // link to dihedral types. All dihedral types should be in use and no
    // non-existent dihedral types should be used.
    std::vector<bool> dihTypesUsedByDataFile;
    dihTypesUsedByDataFile.resize(state_->mdFFDihList_.size(), false);

    fprintf(file, "\r\n");

    for(int i=0; i<(int)dihAtoms1.size(); i++)
    {
        if(dihMDTypeIndices[i] < (int)dihTypesUsedByDataFile.size()) dihTypesUsedByDataFile[dihMDTypeIndices[i]] = true;
        else
        {
            fprintf(file, "\tDihedral type index %i does not exist, there are only %i dihedral types declared in total!\r\n", dihMDTypeIndices[i]+1, (int)dihTypesUsedByDataFile.size());
        }
    }


    // Output information about unused bonds
    for(int i=0; i<(int)bondTypesUsedByDataFile.size(); i++)
    {
        if(!bondTypesUsedByDataFile[i])
        {
            fprintf(file, "\tBond type %i was declared but never used:\r\n", i+1);
            CMDFFBond* bond = state_->mdFFBondList_.get(i);
            if(bond) fprintf(file, "\t\t%s\r\n", bond->getArguments(false).data());
        }
    }


    // Output information about unused anlges
    for(int i=0; i<(int)angleTypesUsedByDataFile.size(); i++)
    {
        if(!angleTypesUsedByDataFile[i])
        {
            fprintf(file, "\tAngle type %i was declared but never used:\r\n", i+1);
            CMDFFAngle* angle = state_->mdFFAngleList_.get(i);
            if(angle) fprintf(file, "\t\t%s\r\n", angle->getArguments(false).data());
        }
    }


    // Output information about unused dihedrals
    for(int i=0; i<(int)dihTypesUsedByDataFile.size(); i++)
    {
        if(!dihTypesUsedByDataFile[i])
        {
            fprintf(file, "\tDihedral type %i was declared but never used:\r\n", i+1);
            CMDFFDih* dih = state_->mdFFDihList_.get(i);
            if(dih) fprintf(file, "\t\t%s\r\n", dih->getArguments(false).data());
        }
    }
}

C3DVector CMolTwisterStateTools::getCenterOfMass(const std::vector<int>& atomIndices, int frameIndex) const
{
    double M = 0.0;
    C3DVector Rc;

    for(int i=0; i<(int)atomIndices.size(); i++)
    {
        double m = state_->atoms_[atomIndices[i]]->m_;
        if(fabs(m) < 1E-5)
        {
            printf("Error: found zero mass!");
            return Rc;
        }
        M+= m;

        if((frameIndex < 0) || (frameIndex >= (int)state_->atoms_[atomIndices[i]]->r_.size()))
        {
            printf("Error: frame %i does not exist!", frameIndex);
            return Rc;
        }

        Rc+= state_->atoms_[atomIndices[i]]->r_[frameIndex]*m;
    }
    if(M != 0.0) Rc*= (1.0 / M);
    else
    {
        printf("Error: found vanishing total mass!");
    }

    return Rc;
}

C3DVector CMolTwisterStateTools::getGeometricCenter(const std::vector<int>& atomIndices, int frameIndex) const
{
    C3DVector Rc;
    int atomCount = (int)atomIndices.size();

    for(int i=0; i<atomCount; i++)
    {
        if((frameIndex < 0) || (frameIndex >= (int)state_->atoms_[atomIndices[i]]->r_.size()))
        {
            printf("Error: frame %i does not exist!", frameIndex);
            return Rc;
        }

        Rc+= state_->atoms_[atomIndices[i]]->r_[frameIndex];
    }
    if(atomCount != 0) Rc*= (1.0 / (double)atomCount);
    else
    {
        printf("Error: found vanishing total atom count!");
    }

    return Rc;
}

C3DVector CMolTwisterStateTools::getGeometricCenter(const std::vector<CAtom*>& atoms, int frameIndex)
{
    C3DVector Rc;
    int atomCount = (int)atoms.size();

    for(int i=0; i<atomCount; i++)
    {
        if((frameIndex < 0) || (frameIndex >= atoms[i]->r_.size()))
        {
            printf("Error: frame %i does not exist!", frameIndex);
            return Rc;
        }

        Rc+= atoms[i]->r_[frameIndex];
    }
    if(atomCount != 0) Rc*= (1.0 / (double)atomCount);
    else
    {
        printf("Error: found vanishing total atom count!");
    }

    return Rc;
}

double CMolTwisterStateTools::measureTotCoulombEnergy(const double* a1to4BondSepCoeffs, int frame) const
{
    double qTot = 0.0;
    double distR12, q1, q2;
    CAtom* atom1Ptr;
    CAtom* atom2Ptr;
    C3DVector R12;
    const double K = 8.98755E9; // 1/(4 PI eps_0) [C^{-2}*N*m^2]
    const double e = 1.602176565E-19; // [C]
    const double NA = 6.02214129E23; // [mol^{-1}]
    const double Ke2 = K * e * e * 1.0E10 / 1000.0 * NA; // [kJAA/mol]

    if(frame == -1) frame = state_->getCurrFrameIndex();

    if(!a1to4BondSepCoeffs)
    {
        for(int i=0; i<(int)state_->atoms_.size(); i++)
        {
            atom1Ptr = state_->atoms_[i].get();
            q1 = atom1Ptr->Q_;
            for(int j=i+1; j<(int)state_->atoms_.size(); j++)
            {
                atom2Ptr = state_->atoms_[j].get();
                q2 = atom2Ptr->Q_;

                R12 = atom2Ptr->r_[frame] - atom1Ptr->r_[frame];
                distR12 = R12.norm();
                if(distR12 == 0.0) distR12 = 1.0E-15;

                qTot+= Ke2*q1*q2 / distR12; // [kJ/mol]
            }
        }
    }
    else
    {
        int bondSep;
        double coeff;

        for(int i=0; i<(int)state_->atoms_.size(); i++)
        {
            atom1Ptr = state_->atoms_[i].get();
            q1 = atom1Ptr->Q_;
            for(int j=i+1; j<(int)state_->atoms_.size(); j++)
            {
                atom2Ptr = state_->atoms_[j].get();
                q2 = atom2Ptr->Q_;

                bondSep = atom1Ptr->getBondSepTo(atom2Ptr);
                if(bondSep == 1) coeff = a1to4BondSepCoeffs[0];
                else if(bondSep == 2)  coeff = a1to4BondSepCoeffs[1];
                else if(bondSep == 3)  coeff = a1to4BondSepCoeffs[2];
                else if(bondSep == 4)  coeff = a1to4BondSepCoeffs[3];
                else                   coeff = 1.0;

                R12 = atom2Ptr->r_[frame] - atom1Ptr->r_[frame];
                distR12 = R12.norm();
                if(distR12 == 0.0) distR12 = 1.0E-15;

                qTot+= coeff * Ke2*q1*q2 / distR12; // [kJ/mol]
            }
        }
    }

    return qTot;
}

double CMolTwisterStateTools::measureCoulombPotential(C3DVector at, const std::vector<std::shared_ptr<CAtom>>* atoms, int frame) const
{
    double qTot = 0.0;
    double distR12, q;
    CAtom* atomPtr;
    C3DVector R12;
    const double K = 8.98755E9; // 1/(4 PI eps_0) [C^{-2}*N*m^2]
    const double e = 1.602176565E-19; // [C]
    const double Ke = K * e * 1.0E10 / 1000.0; // [kJAA/C]

    for(int i=0; i<(int)atoms->size(); i++)
    {
        atomPtr = (*atoms)[i].get();
        q = atomPtr->Q_;

        R12 = atomPtr->r_[frame] - at;
        distR12 = R12.norm();
        if(distR12 == 0.0) distR12 = 1.0E-15;

        qTot+= Ke*q / distR12; // [kJ/C]
    }

    return qTot;
}

double CMolTwisterStateTools::measureCoulombPotential(C3DVector at, int frame) const
{
    if(frame == -1) frame = state_->getCurrFrameIndex();

    return measureCoulombPotential(at, &state_->atoms_, frame);
}

END_CUDA_COMPATIBLE()
