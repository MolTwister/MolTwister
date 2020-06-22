#pragma once
#include "MDFFMatrices.h"
#include "MDForces.h"

class CFunctorCalcForce
{
public:
    CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF);

public:
    void setForceFieldMatrices(const CMDFFMatrices& ffMatrices);
    CMDForces operator()(CMDFFMatrices::CAtom& atom);

private:
    C3DVector calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    C3DVector calcBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    C3DVector calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i);
    C3DVector calcForceBondOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& bondType);
    //    C3DVector calcForceAngle(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j const int& k, const int& i);

private:
    int Natomtypes_;
    int Natoms_;
    int Nbonds_;
    int dim_;
    float Lx_;
    float Ly_;
    float Lz_;
    float cutF_;
    CMDFFMatrices::CAtom* atomList_;
    CMDFFMatrices::CPoint* nonBondFFMatrix_;
    size_t* nonBondFFMatrixFFCount_;
    CMDFFMatrices::CBond* bondList_;
    CMDFFMatrices::CPoint* bondFFList_;
    CMDFFMatrices::CAngle* angleList_;
    CMDFFMatrices::CPoint* angleFFList_;
    CMDFFMatrices::CDihedral* dihedralList_;
    CMDFFMatrices::CPoint* dihedralFFList_;
    CMDFFMatrices::CLastError* lastErrorList_;
};
