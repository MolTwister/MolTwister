#include "MDFFMatrices.h"
#include "FunctorGenCellList.h"
#include "FunctorGenNeighList.h"
#include "../../Tools/MolecularSystemTools.h"
#include "../MDLoop/Printf.h"
#include <functional>

BEGIN_CUDA_COMPATIBLE()

HOST_CALLABLE CMDFFMatrices::CCellList::CCellList()
{
    pbcWidthX_ = 0.0f;
    pbcWidthY_ = 0.0f;
    pbcWidthZ_ = 0.0f;
    pbcLowX_ = 0.0f;
    pbcLowY_ = 0.0f;
    pbcLowZ_ = 0.0f;
    cellCountX_ = 0;
    cellCountY_ = 0;
    cellCountZ_ = 0;
}

void CMDFFMatrices::CCellList::init(CMolTwisterState* state, float rCutoff, float dShell, int numAtoms)
{
    float R = rCutoff + dShell;
    C3DRect pbc = state->view3D_->getPBC();

    pbcWidthX_ = (float)pbc.getWidthX();
    pbcWidthY_ = (float)pbc.getWidthY();
    pbcWidthZ_ = (float)pbc.getWidthZ();

    pbcLowX_ = (float)pbc.rLow_.x_;
    pbcLowY_ = (float)pbc.rLow_.y_;
    pbcLowZ_ = (float)pbc.rLow_.z_;

    cellCountX_ = (int)floor(double(pbcWidthX_ / R));
    cellCountY_ = (int)floor(double(pbcWidthY_ / R));
    cellCountZ_ = (int)floor(double(pbcWidthZ_ / R));

    int totNumCells = cellCountX_ * cellCountY_ * cellCountZ_;

    devCellList_ = mtdevice_vector<int>(numAtoms, -1);
    devCellListEntryPointers_ = mtdevice_vector<CMDFFMatrices::CListPointer>(totNumCells);
    devAtomCellIndices_ = mtdevice_vector<CCellListIndex>(state->atoms_.size());
}

void CMDFFMatrices::CCellList::updatePBC(float Lx, float Ly, float Lz)
{
    pbcWidthX_ = Lx;
    pbcWidthY_ = Ly;
    pbcWidthZ_ = Lz;

    pbcLowX_ = -pbcWidthX_ / 2.0f;
    pbcLowY_ = -pbcWidthY_ / 2.0f;
    pbcLowZ_ = -pbcWidthZ_ / 2.0f;
}

void CMDFFMatrices::CCellList::resetCellList()
{
    for(size_t i=0; i<devCellListEntryPointers_.size(); i++)
    {
        CMDFFMatrices::CListPointer listPointer(-1, 0);
        devCellListEntryPointers_[i] = listPointer;
    }
}

HOST_CALLABLE CMDFFMatrices::CNeighList::CNeighList()
{
    maxNeighbours_ = 0;
    numAtoms_ = 0;
}

void CMDFFMatrices::CNeighList::init(CMolTwisterState* state, int initialMaxNeighbors)
{
    maxNeighbours_ = initialMaxNeighbors;

    size_t numAtoms = state->atoms_.size();
    numAtoms_ = (int)numAtoms;
    devNeighList_ = mtdevice_vector<int>(numAtoms * (size_t)maxNeighbours_);
    devNeighListCount_ = mtdevice_vector<int>(numAtoms);
}

void CMDFFMatrices::CNeighList::resetNeighList()
{
    for(size_t i=0; i<devNeighListCount_.size(); i++)
    {
        devNeighListCount_[i] = 0;
    }
}

void CMDFFMatrices::CNeighList::resizeNeighList(int maxNeighbours)
{
    maxNeighbours_ = maxNeighbours;
    devNeighList_ = mtdevice_vector<int>(numAtoms_ * (size_t)maxNeighbours_);
}

CMDFFMatrices::CMDFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell)
{
    state_ = state;
    rCutoff_ = rCutoff;
    dShell_ = dShell;

    bool bondsAcrossPBC = true;
    prepareFFMatrices(state, stdOut, rCutoff, dShell, bondsAcrossPBC,
                      devAtomList_, devForcesList_,
                      devNonBondFFMatrix_, devNonBondFFMatrixFFCount_,
                      devBondFFList_, devAngleFFList_, devDihedralFFList_,
                      devBondsForAtomLists_, devBondsForAtomListPointers_,
                      devAnglesForAtomLists_, devAnglesForAtomListPointers_,
                      devDihedralsForAtomLists_, devDihedralsForAtomListPointers_,
                      cellList_, neighList_,
                      Natoms_, NatomTypes_, Nbonds_, Nangles_, Ndihedrals_);
}

void CMDFFMatrices::updateAtomList(const mthost_vector<CParticle3D>& atomList)
{
    // Download atom list from host, update it and then upload it to the device
    mthost_vector<CAtom> hostAtomList = devAtomList_;
    if(atomList.size() != hostAtomList.size()) return;

    size_t numAtoms = atomList.size();
    for(size_t i=0; i<numAtoms; i++)
    {
        hostAtomList[i].r_ = atomList[i].r_;
        hostAtomList[i].p_ = atomList[i].p_;
    }

    devAtomList_ = hostAtomList;
}

void CMDFFMatrices::genNeighList(float Lx, float Ly, float Lz)
{
    // Generate cell list
    cellList_.resetCellList();
    cellList_.updatePBC(Lx, Ly, Lz);
    CFunctorGenCellList genCellList;
    genCellList.setForceFieldMatrices(*this);
    mttransform(EXEC_POLICY devAtomList_.begin(), devAtomList_.end(), cellList_.devAtomCellIndices_.begin(), genCellList);
    mtcudaDeviceSynchronize();
    genCellList.assembleCellList(cellList_.devAtomCellIndices_);

    // Generate neighborlists
    int storageRequirement = -1;
    do
    {
        neighList_.resetNeighList();
        CFunctorGenNeighList genNeighList;
        genNeighList.setForceFieldMatrices(*this);
        mttransform(EXEC_POLICY devAtomList_.begin(), devAtomList_.end(), neighList_.devNeighListCount_.begin(), genNeighList);
        mtcudaDeviceSynchronize();
        storageRequirement = genNeighList.checkRequiredMaxNumNeighbours(neighList_.devNeighListCount_);
        if(storageRequirement != -1)
        {
            COut::printf("\tSearching for appropriate neighbour list size [trying %i]\r\n", storageRequirement + 10);
            neighList_.resizeNeighList(storageRequirement + 10);
        }

    } while(storageRequirement != -1);
}

void CMDFFMatrices::prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell, bool bondsAcrossPBC,
                                      mtdevice_vector<CAtom>& devAtomList, mtdevice_vector<CForces>& devForcesList,
                                      mtdevice_vector<CPoint>& devNonBondFFMatrix, mtdevice_vector<size_t>& devNonBondFFMatrixFFCount,
                                      mtdevice_vector<CPoint>& devBondFFList, mtdevice_vector<CPoint>& devAngleFFList, mtdevice_vector<CPoint>& devDihedralFFList,
                                      mtdevice_vector<CBond>& devBondsForAtomLists, mtdevice_vector<CListPointer>& devBondsForAtomListPointers,
                                      mtdevice_vector<CAngle>& devAnglesForAtomLists, mtdevice_vector<CListPointer>& devAnglesForAtomListPointers,
                                      mtdevice_vector<CDihedral>& devDihedralsForAtomLists, mtdevice_vector<CListPointer>& devDihedralsForAtomListPointers,
                                      CCellList& cellList, CNeighList& neighList,
                                      int& Natoms,
                                      int& NatomTypes,
                                      int& Nbonds,
                                      int& Nangles,
                                      int& Ndihedrals) const
{
    const int maxNumFFPerAtomicSet = MAX_FF_PER_ATOMIC_SET;
    const int numPointsInForceProfiles = NUM_POINTS_IN_PROFILES;

    // Generate atom list for GPU and set initial positions
    size_t numAtoms = state->atoms_.size();
    Natoms = (int)numAtoms;
    mthost_vector<CAtom> hostAtomList = mthost_vector<CAtom>(numAtoms);
    mthost_vector<CForces> hostForcesList = mthost_vector<CForces>(numAtoms);
    for(size_t i=0; i<numAtoms; i++)
    {
        if(state->atoms_[i]->r_.size() == 0) continue;
        hostAtomList[i].r_.x_ = state->atoms_[i]->r_[0].x_;
        hostAtomList[i].r_.y_ = state->atoms_[i]->r_[0].y_;
        hostAtomList[i].r_.z_ = state->atoms_[i]->r_[0].z_;
        hostAtomList[i].m_ = (float)state->atoms_[i]->m_;
        hostAtomList[i].q_ = (float)state->atoms_[i]->Q_;
        hostAtomList[i].index_ = (int)i;
    }

    // Obtain bonds specified by MD force field
    std::vector<std::vector<int>> bondDestIndices;
    std::vector<int> molIndices;
    bondDestIndices.resize(state_->atoms_.size());

    std::vector<int> mdBondsFromFF[2];
    std::vector<int> mdTypeIndex;
    CMolTwisterStateTools(state_, stdOut).getAllMDBondsInSystem(mdBondsFromFF[0], mdBondsFromFF[1], mdTypeIndex, bondsAcrossPBC);
    for(int i=0; i<(int)mdBondsFromFF[0].size(); i++)
    {
        if(mdBondsFromFF[0][i] < (int)bondDestIndices.size())
            bondDestIndices[mdBondsFromFF[0][i]].emplace_back(mdBondsFromFF[1][i]);
        if(mdBondsFromFF[1][i] < (int)bondDestIndices.size())
            bondDestIndices[mdBondsFromFF[1][i]].emplace_back(mdBondsFromFF[0][i]);
    }

    CMolecularSystemTools(state_, stdOut).genMolIndices(bondDestIndices, molIndices);
    size_t sizeMolIndicesList = molIndices.size();

    // Assign atom-type indices to each atom in the atom list
    std::vector<std::string> atomTypes;
    state->searchForAtomTypes(atomTypes);
    for(size_t i=0; i<numAtoms; i++)
    {
        hostAtomList[i].typeIndex_ = CMolTwisterState::atomTypeToTypeIndex(atomTypes, state->atoms_[i]->getID());
        if(i < sizeMolIndicesList) hostAtomList[i].molIndex_ = molIndices[i];
        else hostAtomList[i].molIndex_ = -1;
    }

    // Prepare cell list vectors and associated properties
    cellList.init(state, rCutoff, dShell, numAtoms);
    neighList.init(state, 10);

    // Generate non-bonded force-field matrix, [toIndex(row, column, ffIndex, pointIndex)]. Assigned are one ore more 1D force-profiles.
    size_t numAtomTypes = atomTypes.size();
    NatomTypes = (int)numAtomTypes;
    mthost_vector<CPoint> hostNonBondFFMatrix = mthost_vector<CPoint>(numAtomTypes * numAtomTypes * maxNumFFPerAtomicSet * numPointsInForceProfiles);
    mthost_vector<size_t> hostNonBondFFMatrixFFCount = mthost_vector<size_t>(numAtomTypes * numAtomTypes);

    // Set up lambda to convert from pair based plots to CFuncPoint plots
    std::function<std::vector<CPoint>(const std::vector<std::pair<float, float>>&, const std::vector<std::pair<float, float>>&)> toFuncPtPlot
            = [](const std::vector<std::pair<float, float>>& fctInF, const std::vector<std::pair<float, float>>& fctInE)
    {
        std::vector<CPoint> fctOut(fctInF.size());
        for(size_t i=0; i<fctInF.size(); i++) fctOut[i] = CPoint(fctInF[i].first, fctInF[i].second, fctInE[i].second);
        return fctOut;
    };

    // Assign the 1D force-profiles for non-bonded forces
    for(size_t c=0; c<numAtomTypes; c++)
    {
        std::string colTypeName = atomTypes[c];
        for(size_t r=0; r<numAtomTypes; r++)
        {
            std::string rowTypeName = atomTypes[r];

            std::shared_ptr<std::vector<int>> ffIndexList = state->mdFFNonBondedList_.indexFromNames(colTypeName, rowTypeName);
            if(ffIndexList)
            {
                hostNonBondFFMatrixFFCount[toIndexNonBond(r, c, numAtomTypes)] = ffIndexList->size();
                for(size_t i=0; i<ffIndexList->size(); i++)
                {
                    int ffIndex = (*ffIndexList)[i];
                    if(ffIndex >= 0)
                    {
                        CMDFFNonBonded* forceField = state->mdFFNonBondedList_.get(ffIndex);
                        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01f, rCutoff, numPointsInForceProfiles),
                                                                forceField->calc1DPotentialProfile(0.01f, rCutoff, numPointsInForceProfiles));
                        for(size_t j=0; j<plot.size(); j++)
                        {
                            hostNonBondFFMatrix[toIndexNonBond(r, c, ffIndex, j, numAtomTypes, numAtomTypes, maxNumFFPerAtomicSet)] = plot[j];
                        }
                    }
                }
            }
        }
    }

    // Generate bond force field vectors
    std::vector<int> bondAtoms1, bondAtoms2, bondMDTypeIndices;
    CMolTwisterStateTools mtStateTools(state, stdOut);
    mtStateTools.getAllMDBondsInSystem(bondAtoms1, bondAtoms2, bondMDTypeIndices, bondsAcrossPBC);

    mthost_vector<CBond> hostBondList = mthost_vector<CBond>(bondAtoms1.size());
    Nbonds = (int)bondAtoms1.size();
    for(size_t i=0; i<bondAtoms1.size(); i++)
    {
        if(i >= bondAtoms2.size()) continue;
        if(i >= bondMDTypeIndices.size()) continue;

        hostBondList[i].atomIndex1_ = bondAtoms1[i];
        hostBondList[i].atomIndex2_ = bondAtoms2[i];
        hostBondList[i].bondType_ = bondMDTypeIndices[i];
    }

    mthost_vector<CPoint> hostBondFFList = mthost_vector<CPoint>(state->mdFFBondList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFBondList_.size(); i++)
    {
        CMDFFBond* forceField = state->mdFFBondList_.get(i);
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01f, rCutoff, numPointsInForceProfiles),
                                                forceField->calc1DPotentialProfile(0.01f, rCutoff, numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostBondFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate angle force field vectors
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    mtStateTools.getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices, bondsAcrossPBC);

    mthost_vector<CAngle> hostAngleList = mthost_vector<CAngle>(angleAtoms1.size());
    Nangles = (int)angleAtoms1.size();
    for(size_t i=0; i<angleAtoms1.size(); i++)
    {
        if(i >= angleAtoms2.size()) continue;
        if(i >= angleAtoms3.size()) continue;
        if(i >= angleMDTypeIndices.size()) continue;

        hostAngleList[i].atomIndex1_ = angleAtoms1[i];
        hostAngleList[i].atomIndex2_ = angleAtoms2[i];
        hostAngleList[i].atomIndex3_ = angleAtoms3[i];
        hostAngleList[i].angleType_ = angleMDTypeIndices[i];
    }

    mthost_vector<CPoint> hostAngleFFList = mthost_vector<CPoint>(state->mdFFAngleList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFAngleList_.size(); i++)
    {
        CMDFFAngle* forceField = state->mdFFAngleList_.get(i);
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles),
                                                forceField->calc1DPotentialProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostAngleFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate dihedral force field vectors
    std::vector<int> dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices;
    mtStateTools.getAllMDDihedralsInSystem(dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices, bondsAcrossPBC);

    mthost_vector<CDihedral> hostDihedralList = mthost_vector<CDihedral>(dihedralAtoms1.size());
    Ndihedrals = (int)dihedralAtoms1.size();
    for(size_t i=0; i<dihedralAtoms1.size(); i++)
    {
        if(i >= dihedralAtoms2.size()) continue;
        if(i >= dihedralAtoms3.size()) continue;
        if(i >= dihedralAtoms4.size()) continue;
        if(i >= dihedralMDTypeIndices.size()) continue;

        hostDihedralList[i].atomIndex1_ = dihedralAtoms1[i];
        hostDihedralList[i].atomIndex2_ = dihedralAtoms2[i];
        hostDihedralList[i].atomIndex3_ = dihedralAtoms3[i];
        hostDihedralList[i].atomIndex4_ = dihedralAtoms4[i];
        hostDihedralList[i].dihedralType_ = dihedralMDTypeIndices[i];
    }

    mthost_vector<CPoint> hostDihedralFFList = mthost_vector<CPoint>(state->mdFFDihList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFDihList_.size(); i++)
    {
        CMDFFDih* forceField = state->mdFFDihList_.get(i);
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles),
                                                forceField->calc1DPotentialProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostDihedralFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate one list per atom that contains all bonds connected to that atom, where all lists are stored in a single list
    mthost_vector<CBond> hostBondsForAtomLists;
    mthost_vector<CListPointer> hostBondsForAtomListPointers(numAtoms);
    int globIndex = 0;
    for(int k=0; k<(int)numAtoms; k++)
    {
        hostBondsForAtomListPointers[k].indexFirstEntry_ = globIndex;
        int bondCount = 0;
        for(int j=0; j<Nbonds; j++)
        {
            int iBondTo = -1;
            if((int)hostBondList[j].atomIndex1_ == k)
            {
                iBondTo = (int)hostBondList[j].atomIndex2_;
            }
            if((int)hostBondList[j].atomIndex2_ == k)
            {
                iBondTo = (int)hostBondList[j].atomIndex1_;
            }
            if((iBondTo != -1) && (hostBondList[j].bondType_ != -1))
            {
                CBond bond;
                bond.atomIndex1_ = k;
                bond.atomIndex2_ = iBondTo;
                bond.bondType_ = hostBondList[j].bondType_;

                hostBondsForAtomLists.push_back(bond);
                globIndex++;
                bondCount++;
            }
        }
        hostBondsForAtomListPointers[k].numEntries_ = bondCount;
    }

    // Generate one list per atom that contains all angles connected to that atom, where all lists are stored in a single list
    mthost_vector<CAngle> hostAnglesForAtomLists;
    mthost_vector<CListPointer> hostAnglesForAtomListPointers(numAtoms);
    globIndex = 0;
    for(int k=0; k<(int)numAtoms; k++)
    {
        hostAnglesForAtomListPointers[k].indexFirstEntry_ = globIndex;
        int angleCount = 0;
        for(int j=0; j<Nangles; j++)
        {
            int iBondToI = -1;
            int iBondToJ = -1;
            bool kIsCenterOfAngle = false;
            if((int)hostAngleList[j].atomIndex1_ == k)
            {
                iBondToI = (int)hostAngleList[j].atomIndex2_;
                iBondToJ = (int)hostAngleList[j].atomIndex3_;
            }
            if((int)hostAngleList[j].atomIndex2_ == k)
            {
                kIsCenterOfAngle = true;
                iBondToI = (int)hostAngleList[j].atomIndex1_;
                iBondToJ = (int)hostAngleList[j].atomIndex3_;
            }
            if((int)hostAngleList[j].atomIndex3_ == k)
            {
                iBondToI = (int)hostAngleList[j].atomIndex2_;
                iBondToJ = (int)hostAngleList[j].atomIndex1_;
            }
            if((iBondToI != -1) && (iBondToJ != -1) && (hostAngleList[j].angleType_ != -1))
            {
                CAngle angle;
                angle.atomIndex1_ = k;
                angle.atomIndex2_ = iBondToI;
                angle.atomIndex3_ = iBondToJ;
                angle.angleType_ = hostAngleList[j].angleType_;
                angle.assocAtomIsAtCenterOfAngle_ = kIsCenterOfAngle;

                hostAnglesForAtomLists.push_back(angle);
                globIndex++;
                angleCount++;
            }
        }
        hostAnglesForAtomListPointers[k].numEntries_ = angleCount;
    }

    // Generate one list per atom that contains all angles connected to that atom, where all lists are stored in a single list
    mthost_vector<CDihedral> hostDihedralsForAtomLists;
    mthost_vector<CListPointer> hostDihedralsForAtomListPointers(numAtoms);
    globIndex = 0;
    for(int k=0; k<(int)numAtoms; k++)
    {
        hostDihedralsForAtomListPointers[k].indexFirstEntry_ = globIndex;
        int dihedralCount = 0;
        for(int j=0; j<Ndihedrals; j++)
        {
            int iBondToI = -1;
            int iBondToJ = -1;
            int iBondToL = -1;
            bool kIsCenterOfDihedral = false;
            if((int)hostDihedralList[j].atomIndex1_ == k)
            {
                iBondToI = (int)hostDihedralList[j].atomIndex2_;
                iBondToJ = (int)hostDihedralList[j].atomIndex3_;
                iBondToL = (int)hostDihedralList[j].atomIndex4_;
            }
            if((int)hostDihedralList[j].atomIndex2_ == k)
            {
                kIsCenterOfDihedral = true;
                iBondToI = (int)hostDihedralList[j].atomIndex1_;
                iBondToJ = (int)hostDihedralList[j].atomIndex3_;
                iBondToL = (int)hostDihedralList[j].atomIndex4_;
            }
            if((int)hostDihedralList[j].atomIndex3_ == k)
            {
                kIsCenterOfDihedral = true;
                iBondToI = (int)hostDihedralList[j].atomIndex4_;
                iBondToJ = (int)hostDihedralList[j].atomIndex2_;
                iBondToL = (int)hostDihedralList[j].atomIndex1_;
            }
            if((int)hostDihedralList[j].atomIndex4_ == k)
            {
                iBondToI = (int)hostDihedralList[j].atomIndex3_;
                iBondToJ = (int)hostDihedralList[j].atomIndex2_;
                iBondToL = (int)hostDihedralList[j].atomIndex1_;
            }
            if((iBondToI != -1) && (iBondToJ != -1) && (iBondToL != -1) && (hostDihedralList[j].dihedralType_ != -1))
            {
                CDihedral dihedral;
                dihedral.atomIndex1_ = k;
                dihedral.atomIndex2_ = iBondToI;
                dihedral.atomIndex3_ = iBondToJ;
                dihedral.atomIndex4_ = iBondToL;
                dihedral.dihedralType_ = hostDihedralList[j].dihedralType_;
                dihedral.assocAtomIsAtCenterOfDihedral_ = kIsCenterOfDihedral;

                hostDihedralsForAtomLists.push_back(dihedral);
                globIndex++;
                dihedralCount++;
            }
        }
        hostDihedralsForAtomListPointers[k].numEntries_ = dihedralCount;
    }

    // Upload results to device
    devAtomList = hostAtomList;
    devForcesList = hostForcesList;
    devNonBondFFMatrix = hostNonBondFFMatrix;
    devNonBondFFMatrixFFCount = hostNonBondFFMatrixFFCount;
    devBondFFList = hostBondFFList;
    devAngleFFList = hostAngleFFList;
    devDihedralFFList = hostDihedralFFList;
    devBondsForAtomLists = hostBondsForAtomLists;
    devBondsForAtomListPointers = hostBondsForAtomListPointers;
    devAnglesForAtomLists = hostAnglesForAtomLists;
    devAnglesForAtomListPointers = hostAnglesForAtomListPointers;
    devDihedralsForAtomLists = hostDihedralsForAtomLists;
    devDihedralsForAtomListPointers = hostDihedralsForAtomListPointers;
}

HOSTDEV_CALLABLE size_t CMDFFMatrices::toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet)
{
    return columnCount*(rowCount*(maxNumFFPerAtomicSet*pointIndex + ffIndex) + rowIndex) + columnIndex;
}

HOSTDEV_CALLABLE size_t CMDFFMatrices::toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount)
{
    return columnCount*rowIndex + columnIndex;
}

HOSTDEV_CALLABLE size_t CMDFFMatrices::toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles)
{
    return listIndex*numPointsInForceProfiles + pointIndex;
}

END_CUDA_COMPATIBLE()
