#include "FunctorCalcForce.h"
#include <float.h>
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorCalcForce::CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF)
{
    dim_ = dim;
    Lx_ = Lx;
    Ly_ = Ly;
    Lz_ = Lz;
    cutF_ = cutF;

    devAtomList_ = nullptr;
    devNonBondFFMatrix_ = nullptr;
    devNonBondFFMatrixFFCount_ = nullptr;
    devBondList_ = nullptr;
    devBondFFList_ = nullptr;
    devAngleList_ = nullptr;
    devAngleFFList_ = nullptr;
    devDihedralList_ = nullptr;
    devDihedralFFList_ = nullptr;
    devLastErrorList_ = nullptr;
}

void CFunctorCalcForce::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devAtomList_ = mtraw_pointer_cast(&ffMatrices.devAtomList_[0]);
    devNonBondFFMatrix_ = mtraw_pointer_cast(&ffMatrices.devNonBondFFMatrix_[0]);
    devNonBondFFMatrixFFCount_ = mtraw_pointer_cast(&ffMatrices.devNonBondFFMatrixFFCount_[0]);
    devBondList_ = mtraw_pointer_cast(&ffMatrices.devBondList_[0]);
    devBondFFList_ = mtraw_pointer_cast(&ffMatrices.devBondFFList_[0]);
    devAngleList_ = mtraw_pointer_cast(&ffMatrices.devAngleList_[0]);
    devAngleFFList_ = mtraw_pointer_cast(&ffMatrices.devAngleFFList_[0]);
    devDihedralList_ = mtraw_pointer_cast(&ffMatrices.devDihedralList_[0]);
    devDihedralFFList_ = mtraw_pointer_cast(&ffMatrices.devDihedralFFList_[0]);
    devLastErrorList_ = mtraw_pointer_cast(&ffMatrices.devLastErrorList_[0]);

    Natoms_ = ffMatrices.getNumAtoms();
    Natomtypes_ = ffMatrices.getNumAtomTypes();
    Nbonds_ = ffMatrices.getNumBonds();
}

HOSTDEV_CALLABLE CMDFFMatrices::CForces CFunctorCalcForce::operator()(CMDFFMatrices::CAtom& atom)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx_, 0.0,   0.0 );
    C3DVector PBCy = C3DVector( 0.0,  Ly_,  0.0 );
    C3DVector PBCz = C3DVector( 0.0,  0.0,  Lz_ );

    // Clear forces from primary image (pi = primary image)
    C3DVector Fpi(0.0, 0.0, 0.0);

    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    // :TODO: Later this will be a loop over the neighbor list only!!!
    // :TODO: Here, Coulomb is combined into short range. Should make sure that only Coulomb go past PBC!!!
    int k = atom.index_;
    devLastErrorList_[k].reset();
    C3DVector r_k = devAtomList_[k].r_;
    for(int i=0; i<Natoms_; i++)
    {
        C3DVector r_i = devAtomList_[i].r_;

        Fpi+= calcForceNonBondedOn_r_k(r_k, r_i, k, i);

        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCx, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCx, k, i);

        if(dim_ < 2) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCy, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCy, k, i);

        if(dim_ < 3) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCz, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCz, k, i);
    }
    F+= Fpi;

    // Add forces from harmonic bonds on particle k
    for(int j=0; j<Nbonds_; j++)
    {
        int iBondTo = -1;
        if((int)devBondList_[j].atomIndex1_ == k)
        {
            iBondTo = (int)devBondList_[j].atomIndex2_;
        }
        if((int)devBondList_[j].atomIndex2_ == k)
        {
            iBondTo = (int)devBondList_[j].atomIndex1_;
        }
        if(devBondList_[j].bondType_ != -1)
        {
            C3DVector r_k = devAtomList_[k].r_;
            C3DVector r_i = devAtomList_[iBondTo].r_;

            Fpi+= calcForceBondOn_r_k(r_k, r_i, devBondList_[j].bondType_);
        }
    }
    F+= Fpi;

    if(fabs(F.x_) > cutF_) { F.x_ = ((F.x_ >= 0.0) ? 1.0 : -1.0) * cutF_; devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }
    if(fabs(F.y_) > cutF_) { F.y_ = ((F.y_ >= 0.0) ? 1.0 : -1.0) * cutF_; devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }
    if(fabs(F.z_) > cutF_) { F.z_ = ((F.z_ >= 0.0) ? 1.0 : -1.0) * cutF_; devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }

    return CMDFFMatrices::CForces(F, Fpi);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;

    // Return coeffs that correspond to calculating force on r2 (i.e., dr/dx_2, etc.)
    return (r12 * RInv);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;

    // Return coeffs that correspond to calculating force on r2 (i.e., dr/dx_2, etc.)
    return (r12 * RInv);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dr) * (dr/dx_k, dr/dy_k, dr/dz_k)
    // We have (-dU/dr) stored in table form from nonBondFFMatrix_,
    // we just need to retrieve the one stored for (k, i) and interpolate
    // it for r = | r_k - r_i |. We use linear interpolation. If r < 1E-3,
    // we assume that atom i and k are the same in the same image. This,
    // should not yield a contribution to the forces on atom k.
    float r = (r_k - r_i).norm();
    if(r < 1E-3f) return C3DVector(0.0f, 0.0f, 0.0f);

    float mdU_dr_Sum = 0.0f;
    size_t ffCount = devNonBondFFMatrixFFCount_[CMDFFMatrices::toIndexNonBond(devAtomList_[k].typeIndex_, devAtomList_[i].typeIndex_, Natomtypes_)];
    for(int ffIndex=0; ffIndex<(int)ffCount; ffIndex++)
    {
        // Use number of points in plot to get first and last value, thus obtaining r_min and r_max, which
        // is to be used to estimate delta_r, hence enabling calculation of the closest index, j, below r
        CMDFFMatrices::CPoint p_first = devNonBondFFMatrix_[CMDFFMatrices::toIndexNonBond(devAtomList_[k].typeIndex_, devAtomList_[i].typeIndex_, ffIndex, 0, Natomtypes_, Natomtypes_, MAX_FF_PER_ATOMIC_SET)];
        CMDFFMatrices::CPoint p_last = devNonBondFFMatrix_[CMDFFMatrices::toIndexNonBond(devAtomList_[k].typeIndex_, devAtomList_[i].typeIndex_, ffIndex, NUM_POINTS_IN_PROFILES-1, Natomtypes_, Natomtypes_, MAX_FF_PER_ATOMIC_SET)];
        float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
        int idx = int((r - p_first.x_) / delta);

        // Interpolate the function value between i and (i+1)
        // :TODO: We should not include forced beteen same atoms (i.e., if r_k=r_i).
        if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
        {
            CMDFFMatrices::CPoint p1 = devNonBondFFMatrix_[CMDFFMatrices::toIndexNonBond(devAtomList_[k].typeIndex_, devAtomList_[i].typeIndex_, ffIndex, idx, Natomtypes_, Natomtypes_, MAX_FF_PER_ATOMIC_SET)];
            CMDFFMatrices::CPoint p2 = devNonBondFFMatrix_[CMDFFMatrices::toIndexNonBond(devAtomList_[k].typeIndex_, devAtomList_[i].typeIndex_, ffIndex, idx+1, Natomtypes_, Natomtypes_, MAX_FF_PER_ATOMIC_SET)];
            float denom = p1.x_ - p2.x_;
            if(denom == 0.0f) denom = 3.0f*FLT_MIN;
            float a = (p1.y_ - p2.y_) / denom;
            float b = p1.y_ - a*p1.x_;
            mdU_dr_Sum+= (a*r + b);
        }
    }

    // Now that we have (-dU/dr) at r, we find the analytic calculation
    // of dr/dx_k, dr/dy_k and dr/dz_k, where r = sqrt((r_x_k - r_x_i)^2
    // + (r_y_k - r_y_i)^2 + (r_z_k - r_z_i)^2) and then calulcate the
    // forces on r_k.
    C3DVector c = calcNonBondForceCoeffs12(r_i, r_k);
    return c*mdU_dr_Sum;
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceBondOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& bondType)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dr) * (dr/dx_k, dr/dy_k, dr/dz_k)
    // We have (-dU/dr) stored in table form from bondFFList_,
    // we just need to retrieve the one stored for (k, i) and interpolate
    // it for r = | r_k - r_i |. We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining r_min and r_max, which
    // is to be used to estimate delta_r, hence enabling calculation of the closest index, j, below r
    CMDFFMatrices::CPoint p_first = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, 0,NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, NUM_POINTS_IN_PROFILES-1,NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    float r = (r_k - r_i).norm();
    int idx = int((r - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dr = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, idx,NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, idx+1,NUM_POINTS_IN_PROFILES)];
        float denom = p1.x_ - p2.x_;
        if(denom == 0.0f) denom = 3.0f*FLT_MIN;
        float a = (p1.y_ - p2.y_) / denom;
        float b = p1.y_ - a*p1.x_;
        mdU_dr+= (a*r + b);
    }

    // Now that we have (-dU/dr) at r, we find the analytic calculation
    // of dr/dx_k, dr/dy_k and dr/dz_k, where r = sqrt((r_x_k - r_x_i)^2
    // + (r_y_k - r_y_i)^2 + (r_z_k - r_z_i)^2) and then calulcate the
    // forces on r_k.
    C3DVector c = calcBondForceCoeffs12(r_i, r_k);
    return c*mdU_dr;
}

END_CUDA_COMPATIBLE()
