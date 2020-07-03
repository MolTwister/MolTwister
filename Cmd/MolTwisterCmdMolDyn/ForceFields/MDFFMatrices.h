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
    };

    class CCellListIndex
    {
    public:
        HOSTDEV_CALLABLE CCellListIndex() { ix_ = iy_ = iz_ = 0; }
        HOSTDEV_CALLABLE CCellListIndex(int ix, int iy, int iz) { ix_ = ix; iy_ = iy; iz_ = iz; }

    public:
        int ix_;
        int iy_;
        int iz_;
    };

    class CCellList
    {
    public:
        HOST_CALLABLE CCellList();

    public:
        void init(CMolTwisterState* state, float rCutoff, float dShell);
        void resetCellList();

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

    public:
        mtdevice_vector<int> devCellList_;
        mtdevice_vector<int> devCellListCount_;
        mtdevice_vector<CCellListIndex> devAtomCellIndices_;

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
    };

    class CNeighList
    {
    public:
        HOST_CALLABLE CNeighList();

    public:
        void init(CMolTwisterState* state, int maxNeighbors);
        void resetNeighList();
        HOSTDEV_CALLABLE int getMaxNeighbors() const { return maxNeighbours_; }

    public:
        mtdevice_vector<int> devNeighList_;
        mtdevice_vector<int> devNeighListCount_;

    private:
        int maxNeighbours_;
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
    float getRCutoff() const { return rCutoff_; }
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
                           CCellList& cellList, CNeighList& neighList,
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
    CNeighList neighList_;

private:
    float rCutoff_;
    float dShell_;
    int Natoms_;
    int Nbonds_;
    int NatomTypes_;
    CMolTwisterState* state_;
};

END_CUDA_COMPATIBLE()
