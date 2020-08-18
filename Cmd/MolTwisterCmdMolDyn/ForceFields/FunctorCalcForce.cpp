#include "FunctorCalcForce.h"
#include "FunctorGenNeighList.h"
#include <float.h>
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorCalcForce::CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF)
{
    Natomtypes_ = 0;
    Natoms_ = 0;
    Nbonds_ = 0;
    Nangles_ = 0;

    dim_ = dim;
    Lx_ = Lx;
    Ly_ = Ly;
    Lz_ = Lz;
    pbc_ = C3DRect(C3DVector(-(double)Lx/2.0, -(double)Ly/2.0, -(double)Lz/2.0), C3DVector((double)Lx/2.0, (double)Ly/2.0, (double)Lz/2.0));
    cutF_ = cutF;
    maxNeighbours_ = 0;

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
    devNeighList_ = nullptr;
    devNeighListCount_ = nullptr;
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
    devNeighList_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighList_[0]);
    devNeighListCount_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighListCount_[0]);

    Natoms_ = ffMatrices.getNumAtoms();
    Natomtypes_ = ffMatrices.getNumAtomTypes();
    Nbonds_ = ffMatrices.getNumBonds();
    Nangles_ =  ffMatrices.getNumAngles();

    maxNeighbours_ = ffMatrices.neighList_.getMaxNeighbors();
}

HOSTDEV_CALLABLE CMDFFMatrices::CForces CFunctorCalcForce::operator()(CMDFFMatrices::CAtom& atom)
{
    C3DVector F, f;
    C3DVector PBCx = C3DVector( (double)Lx_, 0.0,          0.0         );
    C3DVector PBCy = C3DVector( 0.0,         (double)Ly_,  0.0         );
    C3DVector PBCz = C3DVector( 0.0,         0.0,          (double)Lz_ );

    // Clear forces from primary image (pi = primary image)
    C3DVector Fpi(0.0, 0.0, 0.0);

    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    // :TODO: Later this will be a loop over the neighbor list only!!!
    // :TODO: Here, Coulomb is combined into short range. Should make sure that only Coulomb go past PBC!!!
    int k = atom.index_;
    devLastErrorList_[k].reset();
    C3DVector r_k = devAtomList_[k].r_;
    int numNeighbors = devNeighListCount_[atom.index_];
    for(int neighIndex=0; neighIndex<numNeighbors; neighIndex++)
    {
        int i = devNeighList_[CFunctorGenNeighList::neighIndexToFlatIndex(atom.index_, neighIndex, maxNeighbours_)];
        C3DVector r_i = devAtomList_[i].r_;
        r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);

        f = calcForceNonBondedOn_r_k(r_k, r_i, k, i);
        Fpi+= f;
        F+= f;

        /*
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCx, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCx, k, i);

        if(dim_ < 2) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCy, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCy, k, i);

        if(dim_ < 3) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCz, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCz, k, i);*/
    }

    // Add forces from bonds on particle k
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
        if((iBondTo != -1) && (devBondList_[j].bondType_ != -1))
        {
            C3DVector r_k = devAtomList_[k].r_;
            C3DVector r_i = devAtomList_[iBondTo].r_;
            r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);

            f = calcForceBondOn_r_k(r_k, r_i, devBondList_[j].bondType_);
            Fpi+= f;
            F+= f;
        }
    }

    // Add forces from angles on particle k
    for(int j=0; j<Nangles_; j++)
    {
        int iBondTo2 = -1;
        int iBondTo3 = -1;
        bool kIsCenterOfAngle = false;
        if((int)devAngleList_[j].atomIndex1_ == k)
        {
            iBondTo2 = (int)devAngleList_[j].atomIndex2_;
            iBondTo3 = (int)devAngleList_[j].atomIndex3_;
        }
        if((int)devAngleList_[j].atomIndex2_ == k)
        {
            kIsCenterOfAngle = true;
            iBondTo2 = (int)devAngleList_[j].atomIndex1_;
            iBondTo3 = (int)devAngleList_[j].atomIndex3_;
        }
        if((int)devAngleList_[j].atomIndex3_ == k)
        {
            iBondTo2 = (int)devAngleList_[j].atomIndex1_;
            iBondTo3 = (int)devAngleList_[j].atomIndex2_;
        }
        if((iBondTo2 != -1) && (iBondTo3 != -1) && (devBondList_[j].bondType_ != -1))
        {
            C3DVector r_k = devAtomList_[k].r_;
            C3DVector r_i = devAtomList_[iBondTo2].r_;
            C3DVector r_j = devAtomList_[iBondTo3].r_;
            r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);
            r_k.moveToSameSideOfPBCAsThis(r_j, pbc_);

            f = calcForceAngularOn_r_k(r_k, r_i, devAngleList_[j].angleType_);
            Fpi+= f;
            F+= f;
        }
    }

    if(fabs(F.x_) > double(cutF_)) { F.x_ = ((F.x_ >= 0.0) ? 1.0 : -1.0) * double(cutF_); devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }
    if(fabs(F.y_) > double(cutF_)) { F.y_ = ((F.y_ >= 0.0) ? 1.0 : -1.0) * double(cutF_); devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }
    if(fabs(F.z_) > double(cutF_)) { F.z_ = ((F.z_ >= 0.0) ? 1.0 : -1.0) * double(cutF_); devLastErrorList_[k].lastWarningCode_ = CMDFFMatrices::CLastError::warnForcesWereCut; }

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

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcAngularForceCoeffs13(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3) const
{
    C3DVector   r12 = r1 - r2;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    C3DVector   K = r12*RInv;
    C3DVector   r32 = r3 - r2;
    double      A = K.x_*r32.x_ + K.y_*r32.y_ + K.z_*r32.z_;
    double      B = r32.x_*r32.x_ + r32.y_*r32.y_ + r32.z_*r32.z_;
    double      sqrtB = sqrt(B);
    double      sqrtB3 = (B == 0.0) ? sqrt(1E-30) : sqrtB*sqrtB*sqrtB;
    double      dW_dx = (K.x_*B - r32.x_*A) / sqrtB3;
    double      dW_dy = (K.y_*B - r32.y_*A) / sqrtB3;
    double      dW_dz = (K.z_*B - r32.z_*A) / sqrtB3;
    double      W = (sqrtB == 0.0) ? A / 10E-3 : A / sqrtB;
    double      dDiv = sqrt(1.0 - W*W);
    double      dTheta_dW = (dDiv == 0.0) ? -1.0 / 1E-10 : -1.0 / dDiv;
    C3DVector   dTheta_dr3 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);

    // Return coeffs that correspond to calculating force on r3 (i.e., dtheta/dx_3, etc.)
    return dTheta_dr3;
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcAngularForceCoeffs2(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3) const
{
    C3DVector   r12 = r1 - r2;
    C3DVector   r32 = r3 - r2;
    double      A = r12.x_*r12.x_ + r12.y_*r12.y_ + r12.z_*r12.z_;
    double      B = r32.x_*r32.x_ + r32.y_*r32.y_ + r32.z_*r32.z_;
    double      C = r12.x_*r32.x_ + r12.y_*r32.y_ + r12.z_*r32.z_;
    double      AB = A*B;
    double      sqrtAB = sqrt(AB);
    double      sqrtAB3 = (AB == 0.0) ? sqrt(1E-30) : sqrtAB*sqrtAB*sqrtAB;
    double      dW_dx = ((2.0*r2.x_ - r1.x_ - r3.x_)*AB + (r32.x_*A + r12.x_*B)*C) / sqrtAB3;
    double      dW_dy = ((2.0*r2.y_ - r1.y_ - r3.y_)*AB + (r32.y_*A + r12.y_*B)*C) / sqrtAB3;
    double      dW_dz = ((2.0*r2.z_ - r1.z_ - r3.z_)*AB + (r32.z_*A + r12.z_*B)*C) / sqrtAB3;
    double      W = (sqrtAB == 0.0) ? A / 10E-3 : C / sqrtAB;
    double      dDiv = sqrt(1.0 - W*W);
    double      dTheta_dW = (dDiv == 0.0) ? -1.0 / 1E-10 : -1.0 / dDiv;
    C3DVector   dTheta_dr2 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);

    // Return coeffs that correspond to calculating force on r2 (i.e., dtheta/dx_2, etc.)
    return dTheta_dr2;
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
    float r = (float)(r_k - r_i).norm();
    if(r < 1E-3f) return C3DVector(0.0, 0.0, 0.0);

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
    return c*double(mdU_dr_Sum);
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
    CMDFFMatrices::CPoint p_first = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    float r = (float)(r_k - r_i).norm();
    int idx = int((r - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dr = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, idx, NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devBondFFList_[CMDFFMatrices::toIndexBonded(bondType, idx+1, NUM_POINTS_IN_PROFILES)];
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
    return c*double(mdU_dr);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceAngular13On_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& bondType)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dtheta) * (dtheta/dx_k, dtheta/dy_k, dtheta/dz_k)
    // We have (-dU/dtheta) stored in table form from angleFFList_,
    // we just need to retrieve the one stored for (k, i, j) and interpolate
    // it for theta(r_k, r_i, r_j). We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining r_min and r_max, which
    // is to be used to estimate delta_r, hence enabling calculation of the closest index, j, below r
    CMDFFMatrices::CPoint p_first = devAngleFFList_[CMDFFMatrices::toIndexBonded(bondType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devAngleFFList_[CMDFFMatrices::toIndexBonded(bondType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    // :TODO: Calculate the angle between k, i and j
    float r = (float)(r_k - r_i).norm();
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
    return c*double(mdU_dr);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceAngular2On_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& bondType)
{

}

END_CUDA_COMPATIBLE()
