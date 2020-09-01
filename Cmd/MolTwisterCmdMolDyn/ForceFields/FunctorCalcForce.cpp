#include "FunctorCalcForce.h"
#include "FunctorGenNeighList.h"
#include "../Integrators/Constants.h"
#include <float.h>
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorCalcForce::CFunctorCalcForce(int dim, float Lx, float Ly, float Lz, float cutF, float scale12, float scale13, float scale14)
{
    Natomtypes_ = 0;
    Natoms_ = 0;
    Nbonds_ = 0;
    Nangles_ = 0;
    Ndihedrals_ = 0;

    dim_ = dim;
    Lx_ = Lx;
    Ly_ = Ly;
    Lz_ = Lz;
    pbc_ = C3DRect(C3DVector(-(double)Lx/2.0, -(double)Ly/2.0, -(double)Lz/2.0), C3DVector((double)Lx/2.0, (double)Ly/2.0, (double)Lz/2.0));
    cutF_ = cutF;
    maxNeighbours_ = 0;

    scale12_ = scale12;
    scale13_ = scale13;
    scale14_ = scale14;

    devAtomList_ = nullptr;
    devNonBondFFMatrix_ = nullptr;
    devNonBondFFMatrixFFCount_ = nullptr;
    devBondFFList_ = nullptr;
    devAngleFFList_ = nullptr;
    devDihedralFFList_ = nullptr;
    devLastErrorList_ = nullptr;
    devNeighList_ = nullptr;
    devNeighListCount_ = nullptr;

    devBondsForAtomListPointers_ = nullptr;
    devBondsForAtomLists_ = nullptr;
    devAnglesForAtomListPointers_ = nullptr;
    devAnglesForAtomLists_ = nullptr;
    devDihedralsForAtomListPointers_ = nullptr;
    devDihedralsForAtomLists_ = nullptr;
}

void CFunctorCalcForce::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devAtomList_ = mtraw_pointer_cast(&ffMatrices.devAtomList_[0]);
    devNonBondFFMatrix_ = mtraw_pointer_cast(&ffMatrices.devNonBondFFMatrix_[0]);
    devNonBondFFMatrixFFCount_ = mtraw_pointer_cast(&ffMatrices.devNonBondFFMatrixFFCount_[0]);
    devBondFFList_ = mtraw_pointer_cast(&ffMatrices.devBondFFList_[0]);
    devAngleFFList_ = mtraw_pointer_cast(&ffMatrices.devAngleFFList_[0]);
    devDihedralFFList_ = mtraw_pointer_cast(&ffMatrices.devDihedralFFList_[0]);
    devLastErrorList_ = mtraw_pointer_cast(&ffMatrices.devLastErrorList_[0]);
    devNeighList_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighList_[0]);
    devNeighListCount_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighListCount_[0]);

    devBondsForAtomListPointers_ = mtraw_pointer_cast(&ffMatrices.devBondsForAtomListPointers_[0]);
    devBondsForAtomLists_ = mtraw_pointer_cast(&ffMatrices.devBondsForAtomLists_[0]);
    devAnglesForAtomListPointers_ = mtraw_pointer_cast(&ffMatrices.devAnglesForAtomListPointers_[0]);
    devAnglesForAtomLists_ = mtraw_pointer_cast(&ffMatrices.devAnglesForAtomLists_[0]);
    devDihedralsForAtomListPointers_ = mtraw_pointer_cast(&ffMatrices.devDihedralsForAtomListPointers_[0]);
    devDihedralsForAtomLists_ = mtraw_pointer_cast(&ffMatrices.devDihedralsForAtomLists_[0]);

    Natoms_ = ffMatrices.getNumAtoms();
    Natomtypes_ = ffMatrices.getNumAtomTypes();
    Nbonds_ = ffMatrices.getNumBonds();
    Nangles_ =  ffMatrices.getNumAngles();
    Ndihedrals_ = ffMatrices.getNumDihedrals();

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

    int k = atom.index_;
    devLastErrorList_[k].reset();
    C3DVector r_k = devAtomList_[k].r_;
    float q_k = devAtomList_[k].q_;

    // Add forces from bonds on particle k
    int firstIndexBonds = devBondsForAtomListPointers_[k].indexFirstEntry_;
    int numEntriesBonds = devBondsForAtomListPointers_[k].numEntries_;
    for(int j=0; j<numEntriesBonds; j++)
    {
        CMDFFMatrices::CBond bond = devBondsForAtomLists_[firstIndexBonds + j];
        C3DVector r_i = devAtomList_[bond.atomIndex2_].r_;
        r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);                

        f = calcForceBondOn_r_k(r_k, r_i, bond.bondType_);
        Fpi+= f;
        F+= f;
    }

    // Add forces from angles on particle k
    int firstIndexAngles = devAnglesForAtomListPointers_[k].indexFirstEntry_;
    int numEntriesAngles = devAnglesForAtomListPointers_[k].numEntries_;
    for(int j=0; j<numEntriesAngles; j++)
    {
        CMDFFMatrices::CAngle angle = devAnglesForAtomLists_[firstIndexAngles + j];
        C3DVector r_i = devAtomList_[angle.atomIndex2_].r_;
        C3DVector r_j = devAtomList_[angle.atomIndex3_].r_;
        r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);
        r_k.moveToSameSideOfPBCAsThis(r_j, pbc_);

        if(angle.assocAtomIsAtCenterOfAngle_)
        {
            f = calcForceAngularOn_r_i(r_i, r_k, r_j, angle.angleType_);
        }
        else
        {
            f = calcForceAngularOn_r_k(r_k, r_i, r_j, angle.angleType_);
        }
        Fpi+= f;
        F+= f;
    }

    // Add forces from dihedrals on particle k
    int firstIndexDihedrals = devDihedralsForAtomListPointers_[k].indexFirstEntry_;
    int numEntriesDihedrals = devDihedralsForAtomListPointers_[k].numEntries_;
    for(int j=0; j<numEntriesDihedrals; j++)
    {
        CMDFFMatrices::CDihedral dihedral = devDihedralsForAtomLists_[firstIndexDihedrals + j];
        C3DVector r_i = devAtomList_[dihedral.atomIndex2_].r_;
        C3DVector r_j = devAtomList_[dihedral.atomIndex3_].r_;
        C3DVector r_l = devAtomList_[dihedral.atomIndex4_].r_;
        r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);
        r_k.moveToSameSideOfPBCAsThis(r_j, pbc_);
        r_k.moveToSameSideOfPBCAsThis(r_l, pbc_);

        if(dihedral.assocAtomIsAtCenterOfDihedral_)
        {
            f = calcForceDihedralOn_r_i(r_i, r_k, r_j, r_l, dihedral.dihedralType_);
        }
        else
        {
            f = calcForceDihedralOn_r_k(r_k, r_i, r_j, r_l, dihedral.dihedralType_);
        }
        Fpi+= f;
        F+= f;
    }

    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    // :TODO: Here, Coulomb is combined into short range. Should make sure that only Coulomb go past PBC!!!
    int numNeighbors = devNeighListCount_[atom.index_];
    bool hasBondFactor, hasAngleFactor;
    for(int neighIndex=0; neighIndex<numNeighbors; neighIndex++)
    {
        // Calculate the short-range forces between atom index k and its neareast neighbours, i
        int i = devNeighList_[CFunctorGenNeighList::neighIndexToFlatIndex(atom.index_, neighIndex, maxNeighbours_)];
        C3DVector r_i = devAtomList_[i].r_;
        C3DVector r_i_pbc = r_i;
        r_k.moveToSameSideOfPBCAsThis(r_i, pbc_);
        f = calcForceNonBondedOn_r_k(r_k, r_i, k, i);

        // Calculate the long-range forces between atom index k and its nearest neightbours, i
        float q_i = devAtomList_[i].q_;
        C3DVector r = r_i_pbc - r_k; // :TODO: Check if this should be r_k - r_i!!!
        float R3 = (float)r.norm();
        R3 = R3 * R3 * R3;
        if(R3 != 0.0f)
        {
            f+= ( r * double(( Coulomb_K * q_k * q_i ) / R3) );
        }

        // If the k-i interaction is a 1-2, 1-3 or 1-4 interaction (i.e., one of the intermolecular interactions), then perform an appropriate scaling
        hasBondFactor = false;
        hasAngleFactor = false;

        for(int j=0; j<numEntriesBonds; j++)
        {
            CMDFFMatrices::CBond bond = devBondsForAtomLists_[firstIndexBonds + j];
            if(i == bond.atomIndex2_)
            {
                f*= double(scale12_);
                hasBondFactor = true;
                break;
            }
        }

        if(!hasBondFactor)
        {
            for(int j=0; j<numEntriesAngles; j++)
            {
                CMDFFMatrices::CAngle angle = devAnglesForAtomLists_[firstIndexAngles + j];
                if(!angle.assocAtomIsAtCenterOfAngle_ && (i == angle.atomIndex3_))
                {
                    f*= double(scale13_);
                    hasAngleFactor = true;
                    break;
                }
            }
        }

        if(!hasBondFactor && !hasAngleFactor)
        {
            for(int j=0; j<numEntriesDihedrals; j++)
            {
                CMDFFMatrices::CDihedral dihedral = devDihedralsForAtomLists_[firstIndexDihedrals + j];
                if(!dihedral.assocAtomIsAtCenterOfDihedral_ && (i == dihedral.atomIndex4_))
                {
                    f*= double(scale14_);
                    break;
                }
            }
        }

        // Add the short-range and long-range forces between k and i
        Fpi+= f;
        F+= f;

        /*
        // Apply the minimum image convention for long-range forces
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCx, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCx, k, i);

        if(dim_ < 2) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCy, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCy, k, i);

        if(dim_ < 3) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCz, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCz, k, i);*/
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

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcDihedralForceCoeffs14(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3, const C3DVector& r4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    C3DVector       n1 = r23.cross(r21);
    double          n1Norm = n1.norm();
    double          n1NormInv = (n1Norm == 0.0) ? 1.0 / 1E-10 : 1.0 / n1Norm;
    C3DVector       K = n1 * n1NormInv;
    double          a = r32.y_*r34.z_ - r32.z_*r34.y_;
    double          b = r32.z_*r34.x_ - r32.x_*r34.z_;
    double          c = r32.x_*r34.y_ - r32.y_*r34.x_;
    double          A = K.x_*a + K.y_*b + K.z_*c;
    double          B = a*a + b*b + c*c;
    double          sqrtBInv = (B == 0.0) ? 1.0 / 1E-20 : 1.0 / sqrt(B);
    double          sqrtBInv3 = sqrtBInv*sqrtBInv*sqrtBInv;
    double          W = A*sqrtBInv;
    double          dW_dx = ((K.y_*r32.z_ - K.z_*r32.y_)*B - (r32.z_*(r32.z_*r34.x_ - r32.x_*r34.z_) - r32.y_*(r32.x_*r34.y_ - r32.y_*r34.x_))*A) * sqrtBInv3;
    double          dW_dy = ((K.z_*r32.x_ - K.x_*r32.z_)*B - (r32.x_*(r32.x_*r34.y_ - r32.y_*r34.x_) - r32.z_*(r32.y_*r34.z_ - r32.z_*r34.y_))*A) * sqrtBInv3;
    double          dW_dz = ((K.x_*r32.y_ - K.y_*r32.x_)*B - (r32.y_*(r32.y_*r34.z_ - r32.z_*r34.y_) - r32.x_*(r32.z_*r34.x_ - r32.x_*r34.z_))*A) * sqrtBInv3;
    double          div = sqrt(1.0 - W*W);
    double          dTheta_dW = (div == 0.0) ? 1.0 / 1E-10 : 1.0 / div;
    C3DVector       dTheta_dr4 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);

    // Return coeffs that correspond to calculating force on r4 (i.e., dtheta/dx_4, etc.)
    return dTheta_dr4;
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcDihedralForceCoeffs23(const C3DVector& r1, const C3DVector& r2, const C3DVector& r3, const C3DVector& r4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    double          K1 = r32.z_*r34.x_ - r32.x_*r34.z_;
    double          K2 = r34.z_ - r32.z_;
    double          K3 = r23.z_*r21.x_ - r23.x_*r21.z_;
    double          K4 = r32.x_*r34.y_ - r32.y_*r34.x_;
    double          K5 = r32.y_ - r34.y_;
    double          K6 = r23.x_*r21.y_ - r23.y_*r21.x_;
    double          L1 = r32.y_*r34.z_ - r32.z_*r34.y_;
    double          L2 = -K2;
    double          L3 = r23.y_*r21.z_ - r23.z_*r21.y_;
    double          L4 = K4;
    double          L5 = r34.x_ - r32.x_;
    double          L6 = K6;
    double          M1 = L1;
    double          M2 = -K5;
    double          M3 = L3;
    double          M4 = K1;
    double          M5 = -L5;
    double          M6 = K3;
    double          A = L3*L3 + K3*K3 + K6*K6;
    double          B = L1*L1 + K1*K1 + K4*K4;
    double          C = L3*L1 + K3*K1 + K6*K4;
    double          AB = A*B;
    double          sqrtABInv = (AB == 0.0) ? 1.0 / 1E-20 : 1.0 / sqrt(AB);
    double          sqrtABInv3 = sqrtABInv*sqrtABInv*sqrtABInv;
    double          W = C * sqrtABInv;
    double          dC_dx = K2*K3 + r21.y_*K4 + K5*K6 - r21.z_*K1;
    double          dC_dy = L2*L3 - r21.x_*L4 + L5*L6 + r21.z_*L1;
    double          dC_dz = M2*M3 + r21.x_*M4 + M5*M6 - r21.y_*M1;
    double          dW_dx = (dC_dx*AB - C*((r21.y_*K6 - r21.z_*K3)*B + (K1*K2 + K4*K5)*A)) * sqrtABInv3;
    double          dW_dy = (dC_dy*AB - C*((r21.z_*L3 - r21.x_*L6)*B + (L1*L2 + L4*L5)*A)) * sqrtABInv3;
    double          dW_dz = (dC_dz*AB - C*((r21.x_*M6 - r21.y_*M3)*B + (M1*M2 + M4*M5)*A)) * sqrtABInv3;
    double          div = sqrt(1.0 - W*W);
    double          dTheta_dW = (div == 0.0) ? 1.0 / 1E-10 : 1.0 / div;
    C3DVector       dTheta_dr3 = C3DVector(dTheta_dW*dW_dx, dTheta_dW*dW_dy, dTheta_dW*dW_dz);

    // Return coeffs that correspond to calculating force on r3 (i.e., dtheta/dx_3, etc.)
    return dTheta_dr3;
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

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceAngularOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& angleType)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dtheta) * (dtheta/dx_k, dtheta/dy_k, dtheta/dz_k)
    // We have (-dU/dtheta) stored in table form from angleFFList_,
    // we just need to retrieve the one stored for (k, i, j) and interpolate
    // it for theta(r_k, r_i, r_j). We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining theta_min and theta_max, which
    // is to be used to estimate delta_theta, hence enabling calculation of the closest index, j, below theta
    CMDFFMatrices::CPoint p_first = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    C3DVector r_ik = r_k - r_i;
    C3DVector r_ij = r_j - r_i;
    float dotProd = float(r_ik * r_ij);
    float R_ikR_ij = (float)r_ik.norm() * (float)r_ij.norm();
    float cosAngle = (R_ikR_ij != 0.0f) ? dotProd / R_ikR_ij : 1.0f;
    float theta = (float)acos(double(cosAngle));
    int idx = int((theta - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dtheta = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, idx,NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, idx+1,NUM_POINTS_IN_PROFILES)];
        float denom = p1.x_ - p2.x_;
        if(denom == 0.0f) denom = 3.0f*FLT_MIN;
        float a = (p1.y_ - p2.y_) / denom;
        float b = p1.y_ - a*p1.x_;
        mdU_dtheta+= (a*theta + b);
    }

    // Now that we have (-dU/dtheta) at theta, we find the analytic calculation
    // of dtheta/dx_k, dtheta/dy_k and dtheta/dz_k and then calulcate the
    // forces on r_k.
    C3DVector c = calcAngularForceCoeffs13(r_j, r_i, r_k);
    return c*double(mdU_dtheta);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceAngularOn_r_i(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const int& angleType)
{
    // We know that
    // F_i=-grad_2 U = (-dU/dtheta) * (dtheta/dx_i, dtheta/dy_i, dtheta/dz_i)
    // We have (-dU/dtheta) stored in table form from angleFFList_,
    // we just need to retrieve the one stored for (k, i, j) and interpolate
    // it for theta(r_k, r_i, r_j). We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining theta_min and theta_max, which
    // is to be used to estimate delta_theta, hence enabling calculation of the closest index, j, below theta
    CMDFFMatrices::CPoint p_first = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    C3DVector r_ik = r_k - r_i;
    C3DVector r_ij = r_j - r_i;
    float dotProd = float(r_ik * r_ij);
    float R_ikR_ij = (float)r_ik.norm() * (float)r_ij.norm();
    float cosAngle = (R_ikR_ij != 0.0f) ? dotProd / R_ikR_ij : 1.0f;
    float theta = (float)acos(double(cosAngle));
    int idx = int((theta - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dtheta = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, idx,NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devAngleFFList_[CMDFFMatrices::toIndexBonded(angleType, idx+1,NUM_POINTS_IN_PROFILES)];
        float denom = p1.x_ - p2.x_;
        if(denom == 0.0f) denom = 3.0f*FLT_MIN;
        float a = (p1.y_ - p2.y_) / denom;
        float b = p1.y_ - a*p1.x_;
        mdU_dtheta+= (a*theta + b);
    }

    // Now that we have (-dU/dtheta) at theta, we find the analytic calculation
    // of dtheta/dx_i, dtheta/dy_i and dtheta/dz_i and then calulcate the
    // forces on r_i.
    C3DVector c = calcAngularForceCoeffs2(r_j, r_i, r_k);
    return c*double(mdU_dtheta);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceDihedralOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const C3DVector& r_l, const int& dihedralType)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dtheta) * (dtheta/dx_k, dtheta/dy_k, dtheta/dz_k)
    // We have (-dU/dtheta) stored in table form from dihedralFFList_,
    // we just need to retrieve the one stored for (k, i, j, l) and interpolate
    // it for theta(r_k, r_i, r_j, r_l). We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining theta_min and theta_max, which
    // is to be used to estimate delta_theta, hence enabling calculation of the closest index, j, below theta
    CMDFFMatrices::CPoint p_first = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    C3DVector r_ik = r_k - r_i;
    C3DVector r_ij = r_j - r_i;
    C3DVector r_jl = r_l - r_j;
    C3DVector n1 = r_ik.cross(r_ij);
    C3DVector n2 = r_ij.cross(r_jl);
    float dotProd = float(n1 * n2);
    float n1n2 = (float)n1.norm() * (float)n2.norm();
    float cosAngle = (n1n2 != 0.0f) ? dotProd / n1n2 : 1.0f;
    float theta = (float)acos(double(cosAngle));
    int idx = int((theta - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dtheta = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, idx,NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, idx+1,NUM_POINTS_IN_PROFILES)];
        float denom = p1.x_ - p2.x_;
        if(denom == 0.0f) denom = 3.0f*FLT_MIN;
        float a = (p1.y_ - p2.y_) / denom;
        float b = p1.y_ - a*p1.x_;
        mdU_dtheta+= (a*theta + b);
    }

    // Now that we have (-dU/dtheta) at theta, we find the analytic calculation
    // of dtheta/dx_k, dtheta/dy_k and dtheta/dz_k and then calulcate the
    // forces on r_k.
    C3DVector c = calcDihedralForceCoeffs14(r_l, r_j, r_i, r_k);
    return c*double(mdU_dtheta);
}

HOSTDEV_CALLABLE C3DVector CFunctorCalcForce::calcForceDihedralOn_r_i(const C3DVector& r_k, const C3DVector& r_i, const C3DVector& r_j, const C3DVector& r_l, const int& dihedralType)
{
    // We know that
    // F_i=-grad_1 U = (-dU/dtheta) * (dtheta/dx_i, dtheta/dy_i, dtheta/dz_i)
    // We have (-dU/dtheta) stored in table form from dihedralFFList_,
    // we just need to retrieve the one stored for (k, i, j, l) and interpolate
    // it for theta(r_k, r_i, r_j, r_l). We use linear interpolation.

    // Use number of points in plot to get first and last value, thus obtaining theta_min and theta_max, which
    // is to be used to estimate delta_theta, hence enabling calculation of the closest index, j, below theta
    CMDFFMatrices::CPoint p_first = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, 0, NUM_POINTS_IN_PROFILES)];
    CMDFFMatrices::CPoint p_last = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, NUM_POINTS_IN_PROFILES-1, NUM_POINTS_IN_PROFILES)];
    float delta = (p_last.x_ - p_first.x_) / float(NUM_POINTS_IN_PROFILES);
    C3DVector r_ik = r_k - r_i;
    C3DVector r_ij = r_j - r_i;
    C3DVector r_jl = r_l - r_j;
    C3DVector n1 = r_ik.cross(r_ij);
    C3DVector n2 = r_ij.cross(r_jl);
    float dotProd = float(n1 * n2);
    float n1n2 = (float)n1.norm() * (float)n2.norm();
    float cosAngle = (n1n2 != 0.0f) ? dotProd / n1n2 : 1.0f;
    float theta = (float)acos(double(cosAngle));
    int idx = int((theta - p_first.x_) / delta);

    // Interpolate the function value between i and (i+1)
    float mdU_dtheta = 0.0f;
    if((idx >= 0) && (idx < (NUM_POINTS_IN_PROFILES-1)))
    {
        CMDFFMatrices::CPoint p1 = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, idx,NUM_POINTS_IN_PROFILES)];
        CMDFFMatrices::CPoint p2 = devDihedralFFList_[CMDFFMatrices::toIndexBonded(dihedralType, idx+1,NUM_POINTS_IN_PROFILES)];
        float denom = p1.x_ - p2.x_;
        if(denom == 0.0f) denom = 3.0f*FLT_MIN;
        float a = (p1.y_ - p2.y_) / denom;
        float b = p1.y_ - a*p1.x_;
        mdU_dtheta+= (a*theta + b);
    }

    // Now that we have (-dU/dtheta) at theta, we find the analytic calculation
    // of dtheta/dx_i, dtheta/dy_i and dtheta/dz_i and then calulcate the
    // forces on r_i.
    C3DVector c = calcDihedralForceCoeffs23(r_l, r_j, r_i, r_k);
    return c*double(mdU_dtheta);
}

END_CUDA_COMPATIBLE()
