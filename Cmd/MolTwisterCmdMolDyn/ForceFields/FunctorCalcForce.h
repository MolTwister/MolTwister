#pragma once
#include "MDFFMatrices.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorCalcForce
{
public:
    HOSTDEV_CALLABLE CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF, float cutR, float scale12, float scale13, float scale14, float scale1N);

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE CMDFFMatrices::CForces operator()(CMDFFMatrices::CAtom& atom);

private:
    HOSTDEV_CALLABLE void scaleForcesAndPotentials(C3DVector& f, float& u, const int& atomIndex, const int& molOf_k,
                                                   const int& numEntriesBonds, const int& firstIndexBonds,
                                                   const int& numEntriesAngles, const int& firstIndexAngles,
                                                   const int& numEntriesDihedrals, const int& firstIndexDihedrals);

    HOSTDEV_CALLABLE C3DVector calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    HOSTDEV_CALLABLE C3DVector calcBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    HOSTDEV_CALLABLE C3DVector calcAngularForceCoeffs13(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3) const;
    HOSTDEV_CALLABLE C3DVector calcAngularForceCoeffs2(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3) const;
    HOSTDEV_CALLABLE C3DVector calcDihedralForceCoeffs14(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3, const C3DVector& r4) const;
    HOSTDEV_CALLABLE C3DVector calcDihedralForceCoeffs23(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3, const C3DVector& r4) const;

    HOSTDEV_CALLABLE int calcIndexOfLowerPointOfInterpolationLine(const float& x, const CMDFFMatrices::CPoint& firstPointInList, const CMDFFMatrices::CPoint& lastPointInList) const;
    HOSTDEV_CALLABLE float calcAngle(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j) const;
    HOSTDEV_CALLABLE float calcDihedral(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const C3DVector& r_l) const;
    HOSTDEV_CALLABLE float interpolateLinear(const float& r, const float& x1, const float& y1, const float& x2, const float& y2) const;

    HOSTDEV_CALLABLE C3DVector calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i, float& U) const;
    HOSTDEV_CALLABLE C3DVector calcForceCoulombOn_r_k(const C3DVector& r_k, const C3DVector& r_i, int k, int i, float& U) const;
    HOSTDEV_CALLABLE C3DVector calcForceBondOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& bondType, const int& k, const int& i, float& U) const;
    HOSTDEV_CALLABLE C3DVector calcForceAngularOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& angleType) const;
    HOSTDEV_CALLABLE C3DVector calcForceAngularOn_r_i(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& angleType, float &U) const;
    HOSTDEV_CALLABLE C3DVector calcForceDihedralOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const C3DVector& r_l, const int& dihedralType) const;
    HOSTDEV_CALLABLE C3DVector calcForceDihedralOn_r_i(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const C3DVector& r_l, const int& dihedralType, float &U) const;

private:
    int Natomtypes_;
    int Natoms_;
    int Nbonds_;
    int Nangles_;
    int Ndihedrals_;
    int dim_;
    float Lx_;
    float Ly_;
    float Lz_;
    C3DRect pbc_;
    float cutF_;
    float cutR_;
    float scale12_;
    float scale13_;
    float scale14_;
    float scale1N_;
    CMDFFMatrices::CAtom* devAtomList_;
    CMDFFMatrices::CPoint* devNonBondFFMatrix_;
    size_t* devNonBondFFMatrixFFCount_;
    CMDFFMatrices::CPoint* devBondFFList_;
    CMDFFMatrices::CPoint* devAngleFFList_;
    CMDFFMatrices::CPoint* devDihedralFFList_;
    CMDFFMatrices::CListPointer* devBondsForAtomListPointers_;

    CMDFFMatrices::CBond* devBondsForAtomLists_;
    CMDFFMatrices::CListPointer* devAnglesForAtomListPointers_;
    CMDFFMatrices::CAngle* devAnglesForAtomLists_;
    CMDFFMatrices::CListPointer* devDihedralsForAtomListPointers_;
    CMDFFMatrices::CDihedral* devDihedralsForAtomLists_;

    int* devNeighList_;
    int* devNeighListCount_;
    int maxNeighbours_;
};

END_CUDA_COMPATIBLE()
