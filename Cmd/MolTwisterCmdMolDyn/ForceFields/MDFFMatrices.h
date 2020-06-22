#pragma once
#include "../Integrators/Particle3D.h"
#include "../../Tools/MolTwisterStateTools.h"

#define MAX_FF_PER_ATOMIC_SET 5
#define NUM_POINTS_IN_PROFILES 100

class CMDFFMatrices
{
public:
    class CAtom
    {
    public:
        CAtom() { typeIndex_ = -1; }

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
        CBond() { atomIndex1_ = atomIndex2_ = 0; bondType_ = -1; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        int bondType_;
    };

    class CAngle
    {
    public:
        CAngle() { atomIndex1_ = atomIndex2_ = atomIndex3_ = 0; angleType_ = -1; }

    public:
        size_t atomIndex1_;
        size_t atomIndex2_;
        size_t atomIndex3_;
        int angleType_;
    };

    class CDihedral
    {
    public:
        CDihedral() { atomIndex1_ = atomIndex2_ = atomIndex3_ = atomIndex4_ = 0; dihedralType_ = -1; }

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
        CPoint() { x_ = y_ = 0.0f; }
        CPoint(float x, float y) { x_ = x; y_ = y; }

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
        CLastError() { reset(); }

    public:
        void reset()
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
    void updateAtomList(const std::vector<CParticle3D>& atomList);

    static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet);
    static size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount);
    static size_t toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles);

private:
    void prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, bool bondsAcrossPBC,
                           std::vector<CAtom>& atomList,
                           std::vector<CPoint>& nonBondFFMatrix, std::vector<size_t>& nonBondFFMatrixFFCount,
                           std::vector<CBond>& bondList, std::vector<CPoint>& bondFFList,
                           std::vector<CAngle>& angleList, std::vector<CPoint>& angleFFList,
                           std::vector<CDihedral>& dihedralList, std::vector<CPoint>& dihedralFFList,
                           int& Natoms,
                           int& NatomTypes,
                           int& Nbonds) const;
    void prepareLastErrorList(CMolTwisterState* state, std::vector<CLastError>& lastErrorList) const;

public:
    std::vector<CAtom> atomList_;
    std::vector<CPoint> nonBondFFMatrix_;
    std::vector<size_t> nonBondFFMatrixFFCount_;
    std::vector<CBond> bondList_;
    std::vector<CPoint> bondFFList_;
    std::vector<CAngle> angleList_;
    std::vector<CPoint> angleFFList_;
    std::vector<CDihedral> dihedralList_;
    std::vector<CPoint> dihedralFFList_;
    std::vector<CLastError> lastErrorList_;

private:
    int Natoms_;
    int Nbonds_;
    int NatomTypes_;
    CMolTwisterState* state_;
};
