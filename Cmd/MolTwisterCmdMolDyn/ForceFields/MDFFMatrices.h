#pragma once
#include "../Integrators/Particle3D.h"
#include "../../Tools/MolTwisterStateTools.h"
#include "../../../Utilities/CUDAGeneralizations.h"

#define MAX_FF_PER_ATOMIC_SET 5
#define NUM_POINTS_IN_PROFILES 100

BEGIN_CUDA_COMPATIBLE()

class CMDFFMatrices
{
public:
    class CAtom
    {
    public:
        HOSTDEV_CALLABLE CAtom() { typeIndex_ = -1; }

    public:
        int index_;
        float m_;
        C3DVector r_;
        C3DVector p_;
        int typeIndex_;
        C3DVector* devNeighList_;
        int neighListSize_;
        C3DVector* devNeighListShell_;
        int neighListShellSize_;
    };

    class CCellList
    {
    public:
        HOSTDEV_CALLABLE CCellList()
        {
            pbcWidthX_ = 0;
            pbcWidthY_ = 0;
            pbcWidthZ_ = 0;
            pbcLowX_ = 0;
            pbcLowY_ = 0;
            pbcLowZ_ = 0;
            cellCountX_ = 0;
            cellCountY_ = 0;
            cellCountZ_ = 0;
            maxAtomsInCell_ = 0;
            devCellListRaw_ = nullptr;
            devCellListCountRaw_ = nullptr;
        }

    public:
        void init(CMolTwisterState* state, float rCutoff, float dShell)
        {
            float R = rCutoff + dShell;
            C3DRect pbc = state->view3D_->getPBC();
            maxAtomsInCell_ = ceil(R*R*R);

            pbcWidthX_ = pbc.getWidthX();
            pbcWidthY_ = pbc.getWidthY();
            pbcWidthZ_ = pbc.getWidthZ();

            pbcLowX_ = pbc.rLow_.x_;
            pbcLowY_ = pbc.rLow_.y_;
            pbcLowZ_ = pbc.rLow_.z_;

            cellCountX_ = floor(pbcWidthX_ / R);
            cellCountY_ = floor(pbcWidthY_ / R);
            cellCountZ_ = floor(pbcWidthZ_ / R);

            int totNumCells = cellCountX_ * cellCountY_ * cellCountZ_;

            devCellList_ = mtdevice_vector<int>(totNumCells * maxAtomsInCell_, -1);
            devCellListCount_ = mtdevice_vector<int>(totNumCells, 0);
            devAtomCellIndices_ = mtdevice_vector<int>(state->atoms_.size());

            devCellListRaw_ = mtraw_pointer_cast(&devCellList_[0]);
            devCellListCountRaw_ = mtraw_pointer_cast(&devCellListCount_[0]);
        }

        HOSTDEV_CALLABLE int getPBCWidthX() const { return pbcWidthX_; }
        HOSTDEV_CALLABLE int getPBCWidthY() const { return pbcWidthY_; }
        HOSTDEV_CALLABLE int getPBCWidthZ() const { return pbcWidthZ_; }

        HOSTDEV_CALLABLE int getPBCLowX() const { return pbcLowX_; }
        HOSTDEV_CALLABLE int getPBCLowY() const { return pbcLowY_; }
        HOSTDEV_CALLABLE int getPBCLowZ() const { return pbcLowZ_; }

        HOSTDEV_CALLABLE int getCellCountX() const { return cellCountX_; }
        HOSTDEV_CALLABLE int getCellCountY() const { return cellCountY_; }
        HOSTDEV_CALLABLE int getCellCountZ() const { return cellCountZ_; }

        HOSTDEV_CALLABLE int getMaxAtomsInCell() const { return maxAtomsInCell_; }

        HOSTDEV_CALLABLE int* getCellListRaw() { return devCellListRaw_; }
        HOSTDEV_CALLABLE int* getCellListCountRaw() { return devCellListCountRaw_; }

    public:
        mtdevice_vector<int> devCellList_;
        mtdevice_vector<int> devCellListCount_;
        mtdevice_vector<int> devAtomCellIndices_;

    private:
        int pbcWidthX_;
        int pbcWidthY_;
        int pbcWidthZ_;
        int pbcLowX_;
        int pbcLowY_;
        int pbcLowZ_;
        int cellCountX_;
        int cellCountY_;
        int cellCountZ_;
        int maxAtomsInCell_;
        int* devCellListRaw_;
        int* devCellListCountRaw_;
    };

    class CForces
    {
    public:
        HOSTDEV_CALLABLE CForces() {}
        HOSTDEV_CALLABLE CForces(const C3DVector& F, const C3DVector& Fpi) { F_ = F; Fpi_ = Fpi; }

    public:
        C3DVector F_;
        C3DVector Fpi_;
    };

    class CBond
    {
    public:
        HOSTDEV_CALLABLE CBond() { atomIndex1_ = atomIndex2_ = 0; bondType_ = -1; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        int bondType_;
    };

    class CAngle
    {
    public:
        HOSTDEV_CALLABLE CAngle() { atomIndex1_ = atomIndex2_ = atomIndex3_ = 0; angleType_ = -1; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        size_t atomIndex3_;
        int angleType_;
    };

    class CDihedral
    {
    public:
        HOSTDEV_CALLABLE CDihedral() { atomIndex1_ = atomIndex2_ = atomIndex3_ = atomIndex4_ = 0; dihedralType_ = -1; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        size_t atomIndex3_;
        size_t atomIndex4_;
        int dihedralType_;
    };

    class CPoint
    {
    public:
        HOSTDEV_CALLABLE CPoint() { x_ = y_ = 0.0f; }
        HOSTDEV_CALLABLE CPoint(float x, float y) { x_ = x; y_ = y; }

    public:
        float x_;
        float y_;
    };

    class CLastError
    {
    public:
        enum EErrorCode { errNone = 0 };
        enum EWarningCode { warnNone = 0, warnForcesWereCut };

    public:
        HOSTDEV_CALLABLE CLastError() { reset(); }

    public:
        HOSTDEV_CALLABLE void reset()
        {
            lastErrorCode_ = errNone;
            lastWarningCode_ = warnNone;
        }

    public:
        unsigned char lastErrorCode_;
        unsigned char lastWarningCode_;
    };

public:
    CMDFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell);

public:
    int getNumAtoms() const { return Natoms_; }
    int getNumBonds() const { return Nbonds_; }
    int getNumAtomTypes() const { return NatomTypes_; }
    void updateAtomList(const mthost_vector<CParticle3D>& atomList);
    void genNeighList();

    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet);
    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount);
    HOSTDEV_CALLABLE static size_t toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles);

private:
    void prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell, bool bondsAcrossPBC,
                           mtdevice_vector<CAtom>& devAtomList, mtdevice_vector<CForces>& devForcesList,
                           mtdevice_vector<CPoint>& devNonBondFFMatrix, mtdevice_vector<size_t>& devNonBondFFMatrixFFCount,
                           mtdevice_vector<CBond>& devBondList, mtdevice_vector<CPoint>& devBondFFList,
                           mtdevice_vector<CAngle>& devAngleList, mtdevice_vector<CPoint>& devAngleFFList,
                           mtdevice_vector<CDihedral>& devDihedralList, mtdevice_vector<CPoint>& devDihedralFFList,
                           CCellList& cellList,
                           int& Natoms,
                           int& NatomTypes,
                           int& Nbonds) const;
    void prepareLastErrorList(CMolTwisterState* state, mtdevice_vector<CLastError>& devLastErrorList) const;

public:
    mtdevice_vector<CAtom> devAtomList_;
    mtdevice_vector<CForces> devForcesList_;
    mtdevice_vector<CPoint> devNonBondFFMatrix_;
    mtdevice_vector<size_t> devNonBondFFMatrixFFCount_;
    mtdevice_vector<CBond> devBondList_;
    mtdevice_vector<CPoint> devBondFFList_;
    mtdevice_vector<CAngle> devAngleList_;
    mtdevice_vector<CPoint> devAngleFFList_;
    mtdevice_vector<CDihedral> devDihedralList_;
    mtdevice_vector<CPoint> devDihedralFFList_;
    mtdevice_vector<CLastError> devLastErrorList_;
    CCellList cellList_;

private:
    float rCutoff_;
    float dShell_;
    int Natoms_;
    int Nbonds_;
    int NatomTypes_;
    CMolTwisterState* state_;
};

END_CUDA_COMPATIBLE()
