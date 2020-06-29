#include "MDFFMatrices.h"
#include <functional>

BEGIN_CUDA_COMPATIBLE()

CMDFFMatrices::CMDFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff)
{
    state_ = state;

    bool bondsAcrossPBC = true;
    prepareFFMatrices(state, stdOut, rCutoff, bondsAcrossPBC,
                      devAtomList_, devForcesList_,
                      devNonBondFFMatrix_, devNonBondFFMatrixFFCount_,
                      devBondList_, devBondFFList_,
                      devAngleList_, devAngleFFList_,
                      devDihedralList_, devDihedralFFList_,
                      Natoms_, NatomTypes_, Nbonds_);
    prepareLastErrorList(state, devLastErrorList_);
}

void CMDFFMatrices::updateAtomList(const mthost_vector<CParticle3D>& atomList)
{
    // Download atom list from host, update it and then upload it to the device
    mthost_vector<CAtom> hostAtomList = devAtomList_;
    if(atomList.size() != hostAtomList.size()) return;

    size_t numAtoms = atomList.size();
    for(size_t i=0; i<numAtoms; i++)
    {
        hostAtomList[i].r_ = atomList[i].x;
        hostAtomList[i].p_ = atomList[i].p;
    }

    devAtomList_ = hostAtomList;
}

void CMDFFMatrices::prepareLastErrorList(CMolTwisterState* state, mtdevice_vector<CLastError>& devLastErrorList) const
{
    size_t numAtoms = state->atoms_.size();
    mthost_vector<CLastError> hostLastErrorList = mthost_vector<CLastError>(numAtoms);

    devLastErrorList = hostLastErrorList;
}

void CMDFFMatrices::prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, bool bondsAcrossPBC,
                                      mtdevice_vector<CAtom>& devAtomList, mtdevice_vector<CForces>& devForcesList,
                                      mtdevice_vector<CPoint>& devNonBondFFMatrix, mtdevice_vector<size_t>& devNonBondFFMatrixFFCount,
                                      mtdevice_vector<CBond>& devBondList, mtdevice_vector<CPoint>& devBondFFList,
                                      mtdevice_vector<CAngle>& devAngleList, mtdevice_vector<CPoint>& devAngleFFList,
                                      mtdevice_vector<CDihedral>& devDihedralList, mtdevice_vector<CPoint>& devDihedralFFList,
                                      int& Natoms,
                                      int& NatomTypes,
                                      int& Nbonds) const
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
        hostAtomList[i].m_ = state->atoms_[i]->m_;
        hostAtomList[i].index_ = (int)i;
    }

    // Assign atom-type indices to each atom in the atom list
    std::vector<std::string> atomTypes;
    state->searchForAtomTypes(atomTypes);
    for(size_t i=0; i<numAtoms; i++)
    {
        hostAtomList[i].typeIndex_ = CMolTwisterState::atomTypeToTypeIndex(atomTypes, state->atoms_[i]->getID());
    }

    // Generate non-bonded force-field matrix, [toIndex(row, column, ffIndex, pointIndex)]. Assigned are one ore more 1D force-profiles.
    size_t numAtomTypes = atomTypes.size();
    NatomTypes = (int)numAtomTypes;
    mthost_vector<CPoint> hostNonBondFFMatrix = mthost_vector<CPoint>(numAtomTypes * numAtomTypes * maxNumFFPerAtomicSet * numPointsInForceProfiles);
    mthost_vector<size_t> hostNonBondFFMatrixFFCount = mthost_vector<size_t>(numAtomTypes * numAtomTypes);

    // Set up lambda to convert from pair based plots to CFuncPoint plots
    std::function<std::vector<CPoint>(const std::vector<std::pair<float, float>>&)> toFuncPtPlot = [](const std::vector<std::pair<float, float>>& fctIn)
    {
        std::vector<CPoint> fctOut(fctIn.size());
        for(size_t i=0; i<fctIn.size(); i++) fctOut[i] = CPoint(fctIn[i].first, fctIn[i].second);
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
                        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01, rCutoff, numPointsInForceProfiles));
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
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01, rCutoff, numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostBondFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate angle force field vectors
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    mtStateTools.getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices, bondsAcrossPBC);

    mthost_vector<CAngle> hostAngleList = mthost_vector<CAngle>(angleAtoms1.size());
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
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostAngleFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate dihedral force field vectors
    std::vector<int> dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices;
    mtStateTools.getAllMDDihedralsInSystem(dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices, bondsAcrossPBC);

    mthost_vector<CDihedral> hostDihedralList = mthost_vector<CDihedral>(dihedralAtoms1.size());
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
        std::vector<CPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            hostDihedralFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Upload results to device
    devAtomList = hostAtomList;
    devForcesList = hostForcesList;
    devNonBondFFMatrix = hostNonBondFFMatrix;
    devNonBondFFMatrixFFCount = hostNonBondFFMatrixFFCount;
    devBondList = hostBondList;
    devBondFFList = hostBondFFList;
    devAngleList = hostAngleList;
    devAngleFFList = hostAngleFFList;
    devDihedralList = hostDihedralList;
    devDihedralFFList = hostDihedralFFList;
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
