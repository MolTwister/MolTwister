#pragma once
#include "MDFFMatrices.h"
#include "MDForces.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorCalcForce
{
public:
    HOSTDEV_CALLABLE CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF);

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE CMDForces operator()(CMDFFMatrices::CAtom& atom);

private:
    HOSTDEV_CALLABLE C3DVector calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    HOSTDEV_CALLABLE C3DVector calcBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    HOSTDEV_CALLABLE C3DVector calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i);
    HOSTDEV_CALLABLE C3DVector calcForceBondOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& bondType);
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
    CMDFFMatrices::CAtom* devAtomList_;
    CMDFFMatrices::CPoint* devNonBondFFMatrix_;
    size_t* devNonBondFFMatrixFFCount_;
    CMDFFMatrices::CBond* devBondList_;
    CMDFFMatrices::CPoint* devBondFFList_;
    CMDFFMatrices::CAngle* devAngleList_;
    CMDFFMatrices::CPoint* devAngleFFList_;
    CMDFFMatrices::CDihedral* devDihedralList_;
    CMDFFMatrices::CPoint* devDihedralFFList_;
    CMDFFMatrices::CLastError* devLastErrorList_;
};

END_CUDA_COMPATIBLE()
