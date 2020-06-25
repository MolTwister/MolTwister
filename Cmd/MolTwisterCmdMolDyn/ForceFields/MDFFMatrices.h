#pragma once
#include "../Integrators/Particle3D.h"
#include "../../Tools/MolTwisterStateTools.h"
#include "../Common/AlgorithmDefs.h"

#define MAX_FF_PER_ATOMIC_SET 5
#define NUM_POINTS_IN_PROFILES 100

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
        C3DVector F_;
        C3DVector Fpi_;
        // :TODO: Neighboorlist
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
    CMDFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff);

public:
    int getNumAtoms() const { return Natoms_; }
    int getNumBonds() const { return Nbonds_; }
    int getNumAtomTypes() const { return NatomTypes_; }
    void updateAtomList(const mthost_vector<CParticle3D>& atomList);

    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet);
    HOSTDEV_CALLABLE static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount);
    HOSTDEV_CALLABLE static size_t toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles);

private:
    void prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, bool bondsAcrossPBC,
                           mtdevice_vector<CAtom>& atomList,
                           mtdevice_vector<CPoint>& nonBondFFMatrix, mtdevice_vector<size_t>& nonBondFFMatrixFFCount,
                           mtdevice_vector<CBond>& bondList, mtdevice_vector<CPoint>& bondFFList,
                           mtdevice_vector<CAngle>& angleList, mtdevice_vector<CPoint>& angleFFList,
                           mtdevice_vector<CDihedral>& dihedralList, mtdevice_vector<CPoint>& dihedralFFList,
                           int& Natoms,
                           int& NatomTypes,
                           int& Nbonds) const;
    void prepareLastErrorList(CMolTwisterState* state, mtdevice_vector<CLastError>& lastErrorList) const;

public:
    mtdevice_vector<CAtom> atomList_;
    mtdevice_vector<CPoint> nonBondFFMatrix_;
    mtdevice_vector<size_t> nonBondFFMatrixFFCount_;
    mtdevice_vector<CBond> bondList_;
    mtdevice_vector<CPoint> bondFFList_;
    mtdevice_vector<CAngle> angleList_;
    mtdevice_vector<CPoint> angleFFList_;
    mtdevice_vector<CDihedral> dihedralList_;
    mtdevice_vector<CPoint> dihedralFFList_;
    mtdevice_vector<CLastError> lastErrorList_;

private:
    int Natoms_;
    int Nbonds_;
    int NatomTypes_;
    CMolTwisterState* state_;
};
