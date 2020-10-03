#pragma once
#include "../Integrators/Particle3D.h"
#include "../../Tools/MolTwisterStateTools.h"
#include "../../../Utilities/CUDAGeneralizations.h"

#define MAX_FF_PER_ATOMIC_SET 5
#define NUM_POINTS_IN_PROFILES 1000

BEGIN_CUDA_COMPATIBLE()

class CMDFFMatrices
{
public:
    class CListPointer
    {
    public:
        HOSTDEV_CALLABLE CListPointer() { indexFirstEntry_ = -1; numEntries_ = 0; }
        HOSTDEV_CALLABLE CListPointer(int indexFirstEntry, int numEntries) { indexFirstEntry_ = indexFirstEntry; numEntries_ = numEntries; }

    public:
        int numEntries_;
        int indexFirstEntry_;
    };

    class CAtom
    {
    public:
        HOSTDEV_CALLABLE CAtom() { typeIndex_ = -1; }

    public:
        int index_;
        float m_;
        float q_;
        C3DVector r_;
        C3DVector p_;
        int typeIndex_;
        int molIndex_;
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
        void init(CMolTwisterState* state, float rCutoff, float dShell, int numAtoms);
        void resetCellList();
        void updatePBC(float Lx, float Ly, float Lz);

        HOSTDEV_CALLABLE float getPBCWidthX() const { return pbcWidthX_; }
        HOSTDEV_CALLABLE float getPBCWidthY() const { return pbcWidthY_; }
        HOSTDEV_CALLABLE float getPBCWidthZ() const { return pbcWidthZ_; }

        HOSTDEV_CALLABLE float getPBCLowX() const { return pbcLowX_; }
        HOSTDEV_CALLABLE float getPBCLowY() const { return pbcLowY_; }
        HOSTDEV_CALLABLE float getPBCLowZ() const { return pbcLowZ_; }

        HOSTDEV_CALLABLE int getCellCountX() const { return cellCountX_; }
        HOSTDEV_CALLABLE int getCellCountY() const { return cellCountY_; }
        HOSTDEV_CALLABLE int getCellCountZ() const { return cellCountZ_; }

    public:
        mtdevice_vector<int> devCellList_;
        mtdevice_vector<CListPointer> devCellListEntryPointers_;
        mtdevice_vector<CCellListIndex> devAtomCellIndices_;

    private:
        float pbcWidthX_;
        float pbcWidthY_;
        float pbcWidthZ_;
        float pbcLowX_;
        float pbcLowY_;
        float pbcLowZ_;
        int cellCountX_;
        int cellCountY_;
        int cellCountZ_;
    };

    class CNeighList
    {
    public:
        HOST_CALLABLE CNeighList();

    public:
        void init(CMolTwisterState* state, int initialMaxNeighbors);
        void resetNeighList();
        void resizeNeighList(int maxNeighbours);
        HOSTDEV_CALLABLE int getMaxNeighbors() const { return maxNeighbours_; }

    public:
        mtdevice_vector<int> devNeighList_;
        mtdevice_vector<int> devNeighListCount_;

    private:
        int maxNeighbours_;
        int numAtoms_;
    };

    class CForces
    {
    public:
        HOSTDEV_CALLABLE CForces() { U_ = 0.0f; }
        HOSTDEV_CALLABLE CForces(const C3DVector& F, const C3DVector& Fpi) { F_ = F; Fpi_ = Fpi; U_ = 0.0f; }
        HOSTDEV_CALLABLE CForces(const C3DVector& F, const C3DVector& Fpi, const float& U) { F_ = F; Fpi_ = Fpi; U_ = U; }

    public:
        C3DVector F_;
        C3DVector Fpi_;
        float U_;
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
        HOSTDEV_CALLABLE CAngle() { atomIndex1_ = atomIndex2_ = atomIndex3_ = 0; angleType_ = -1; assocAtomIsAtCenterOfAngle_ = false; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        size_t atomIndex3_;
        int angleType_;
        bool assocAtomIsAtCenterOfAngle_;
    };

    class CDihedral
    {
    public:
        HOSTDEV_CALLABLE CDihedral() { atomIndex1_ = atomIndex2_ = atomIndex3_ = atomIndex4_ = 0; dihedralType_ = -1; assocAtomIsAtCenterOfDihedral_ = false; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        size_t atomIndex3_;
        size_t atomIndex4_;
        int dihedralType_;
        bool assocAtomIsAtCenterOfDihedral_;
    };

    class CPoint
    {
    public:
        HOSTDEV_CALLABLE CPoint() { x_ = f_ = e_ = 0.0f; }
        HOSTDEV_CALLABLE CPoint(float x, float f, float e) { x_ = x; f_ = f; e_ = e; }

    public:
        float x_;
        float f_;
        float e_;
    };

public:
    CMDFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell);

public:
    int getNumAtoms() const { return Natoms_; }
    int getNumBonds() const { return Nbonds_; }
    int getNumAngles() const { return Nangles_; }
    int getNumDihedrals() const { return Ndihedrals_; }
    int getNumAtomTypes() const { return NatomTypes_; }
    float getRCutoff() const { return rCutoff_; }
    void updateAtomList(const mthost_vector<CParticle3D>& atomList);
    void genNeighList(float Lx, float Ly, float Lz);

    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet);
    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount);
    HOSTDEV_CALLABLE static size_t toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles);

private:
    void prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, float dShell, bool bondsAcrossPBC,
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
                           int& Ndihedrals) const;

public:
    mtdevice_vector<CAtom> devAtomList_;
    mtdevice_vector<CForces> devForcesList_;
    mtdevice_vector<CPoint> devNonBondFFMatrix_;
    mtdevice_vector<size_t> devNonBondFFMatrixFFCount_;
    mtdevice_vector<CPoint> devBondFFList_;
    mtdevice_vector<CPoint> devAngleFFList_;
    mtdevice_vector<CPoint> devDihedralFFList_;
    mtdevice_vector<CListPointer> devBondsForAtomListPointers_;
    mtdevice_vector<CBond> devBondsForAtomLists_;
    mtdevice_vector<CListPointer> devAnglesForAtomListPointers_;
    mtdevice_vector<CAngle> devAnglesForAtomLists_;
    mtdevice_vector<CListPointer> devDihedralsForAtomListPointers_;
    mtdevice_vector<CDihedral> devDihedralsForAtomLists_;
    CCellList cellList_;
    CNeighList neighList_;

private:
    float rCutoff_;
    float dShell_;
    int Natoms_;
    int Nbonds_;
    int Nangles_;
    int Ndihedrals_;
    int NatomTypes_;
    CMolTwisterState* state_;
};

END_CUDA_COMPATIBLE()
