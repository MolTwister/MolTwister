#include "VelVerlet.h"
#include "../MDLoop/Printf.h"
#include "Math.h"
#include <float.h>
#include <math.h>
#include <functional>
#include "../../Tools/MolTwisterStateTools.h"

#define MAX_FF_PER_ATOMIC_SET 5

// :TODO: Move all the below class definitions, and their implementations, into its own *.h/*.cu module(s) once they become CUDA compatible
template<class T> T* raw_pointer_cast(const T* ptr) { return (T*)ptr; } // :TODO: Do I need this one?

size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet)
{
    return columnCount*(rowCount*(maxNumFFPerAtomicSet*pointIndex + ffIndex) + rowIndex) + columnIndex;
}

size_t toIndexNonBond(size_t rowIndex, size_t columnIndex, size_t columnCount)
{
    return columnCount*rowIndex + columnIndex;
}

size_t toIndexBonded(size_t listIndex, size_t pointIndex, size_t numPointsInForceProfiles)
{
    return listIndex*numPointsInForceProfiles + pointIndex;
}

class CDevAtom
{
public:
    CDevAtom() { typeIndex_ = -1; }

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

class CDevBond
{
public:
    CDevBond() { atomIndex1_ = atomIndex2_ = 0; bondType_ = -1; }

public:
    size_t atomIndex1_;
    size_t atomIndex2_;
    int bondType_;
};

class CDevAngle
{
public:
    CDevAngle() { atomIndex1_ = atomIndex2_ = atomIndex3_ = 0; angleType_ = -1; }

public:
    size_t atomIndex1_;
    size_t atomIndex2_;
    size_t atomIndex3_;
    int angleType_;
};

class CDevDihedral
{
public:
    CDevDihedral() { atomIndex1_ = atomIndex2_ = atomIndex3_ = atomIndex4_ = 0; dihedralType_ = -1; }

public:
    size_t atomIndex1_;
    size_t atomIndex2_;
    size_t atomIndex3_;
    size_t atomIndex4_;
    int dihedralType_;
};

class CFunctPoint
{
public:
    CFunctPoint() { x_ = y_ = 0.0f; }
    CFunctPoint(float x, float y) { x_ = x; y = y_; }

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

class CDevForceFieldMatrices
{
public:
    CDevForceFieldMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff);

private:
    void prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, bool bondsAcrossPBC,
                           std::vector<CDevAtom>& atomList,
                           std::vector<CFunctPoint>& nonBondFFMatrix, std::vector<size_t>& nonBondFFMatrixFFCount,
                           std::vector<CDevBond>& bondList, std::vector<CFunctPoint>& bondFFList,
                           std::vector<CDevAngle>& angleList, std::vector<CFunctPoint>& angleFFList,
                           std::vector<CDevDihedral>& dihedralList, std::vector<CFunctPoint>& dihedralFFList) const;
    void prepareLastErrorList(CMolTwisterState* state, std::vector<CLastError>& lastErrorList) const;

public:
    std::vector<CDevAtom> atomList_;
    std::vector<CFunctPoint> nonBondFFMatrix_;
    std::vector<size_t> nonBondFFMatrixFFCount_;
    std::vector<CDevBond> bondList_;
    std::vector<CFunctPoint> bondFFList_;
    std::vector<CDevAngle> angleList_;
    std::vector<CFunctPoint> angleFFList_;
    std::vector<CDevDihedral> dihedralList_;
    std::vector<CFunctPoint> dihedralFFList_;
    std::vector<CLastError> lastErrorList_;
};

CDevForceFieldMatrices::CDevForceFieldMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff)
{
    bool bondsAcrossPBC = true;
    prepareFFMatrices(state, stdOut, rCutoff, bondsAcrossPBC, atomList_,
                      nonBondFFMatrix_, nonBondFFMatrixFFCount_,
                      bondList_, bondFFList_,
                      angleList_, angleFFList_,
                      dihedralList_, dihedralFFList_);
    prepareLastErrorList(state, lastErrorList_);
}

void CDevForceFieldMatrices::prepareLastErrorList(CMolTwisterState* state, std::vector<CLastError>& lastErrorList) const
{
    size_t numAtoms = state->atoms_.size();
    lastErrorList = std::vector<CLastError>(numAtoms);
}

void CDevForceFieldMatrices::prepareFFMatrices(CMolTwisterState* state, FILE* stdOut, float rCutoff, bool bondsAcrossPBC,
                                      std::vector<CDevAtom>& atomList,
                                      std::vector<CFunctPoint>& nonBondFFMatrix, std::vector<size_t>& nonBondFFMatrixFFCount,
                                      std::vector<CDevBond>& bondList, std::vector<CFunctPoint>& bondFFList,
                                      std::vector<CDevAngle>& angleList, std::vector<CFunctPoint>& angleFFList,
                                      std::vector<CDevDihedral>& dihedralList, std::vector<CFunctPoint>& dihedralFFList) const
{
    const int maxNumFFPerAtomicSet = MAX_FF_PER_ATOMIC_SET;
    const int numPointsInForceProfiles = 100;

    // Generate atom list for GPU and set initial positions
    size_t numAtoms = state->atoms_.size();
    atomList = std::vector<CDevAtom>(numAtoms);
    for(size_t i=0; i<numAtoms; i++)
    {
        if(state->atoms_[i]->r_.size() == 0) continue;
        atomList[i].r_.x_ = state->atoms_[i]->r_[0].x_;
        atomList[i].r_.y_ = state->atoms_[i]->r_[0].y_;
        atomList[i].r_.z_ = state->atoms_[i]->r_[0].z_;
        atomList[i].m_ = state->atoms_[i]->m_;
        atomList[i].index_ = (int)i;
    }

    // Assign atom-type indices to each atom in the atom list
    std::vector<std::string> atomTypes;
    state->searchForAtomTypes(atomTypes);
    for(size_t i=0; i<numAtoms; i++)
    {
        atomList[i].typeIndex_ = CMolTwisterState::atomTypeToTypeIndex(atomTypes, state->atoms_[i]->getID());
    }

    // Generate non-bonded force-field matrix, [toIndex(row, column, ffIndex, pointIndex)]. Assigned are one ore more 1D force-profiles.
    size_t numAtomTypes = atomTypes.size();
    nonBondFFMatrix = std::vector<CFunctPoint>(numAtomTypes * numAtomTypes * maxNumFFPerAtomicSet * numPointsInForceProfiles);
    nonBondFFMatrixFFCount = std::vector<size_t>(numAtomTypes * numAtomTypes);

    // Set up lambda to convert from pair based plots to CFuncPoint plots
    std::function<std::vector<CFunctPoint>(const std::vector<std::pair<float, float>>&)> toFuncPtPlot = [](const std::vector<std::pair<float, float>>& fctIn)
    {
        std::vector<CFunctPoint> fctOut(fctIn.size());
        for(size_t i=0; i<fctIn.size(); i++) fctOut[i] = CFunctPoint(fctIn[i].first, fctIn[i].second);
        return fctOut;
    };

    // Assign the 1D force-profiles for non-bonded forces
    for(size_t c=0; c<numAtomTypes; c++)
    {
        std::string colTypeName = atomTypes[c];
        for(size_t r=0; r<numAtomTypes; r++)
        {
            std::string rowTypeName = atomTypes[r];

            std::shared_ptr<std::vector<int>> ffIndexList = state->mdFFNonBondedList_.indexFromNames(colTypeName, rowTypeName);
            if(ffIndexList)
            {
                nonBondFFMatrixFFCount[toIndexNonBond(r, c, numAtomTypes)] = ffIndexList->size();
                for(size_t i=0; i<ffIndexList->size(); i++)
                {
                    int ffIndex = (*ffIndexList)[i];
                    if(ffIndex >= 0)
                    {
                        CMDFFNonBonded* forceField = state->mdFFNonBondedList_.get(ffIndex);
                        std::vector<CFunctPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01, rCutoff, numPointsInForceProfiles));
                        for(size_t j=0; j<plot.size(); j++)
                        {
                            nonBondFFMatrix[toIndexNonBond(r, c, ffIndex, j, numAtomTypes, numAtomTypes, maxNumFFPerAtomicSet)] = plot[j];
                        }
                    }
                }
            }
        }
    }

    // Generate bond force field vectors
    std::vector<int> bondAtoms1, bondAtoms2, bondMDTypeIndices;
    CMolTwisterStateTools mtStateTools(state, stdOut);
    mtStateTools.getAllMDBondsInSystem(bondAtoms1, bondAtoms2, bondMDTypeIndices, bondsAcrossPBC);

    bondList = std::vector<CDevBond>(bondAtoms1.size());
    for(size_t i=0; i<bondAtoms1.size(); i++)
    {
        if(i >= bondAtoms2.size()) continue;
        if(i >= bondMDTypeIndices.size()) continue;

        bondList[i].atomIndex1_ = bondAtoms1[i];
        bondList[i].atomIndex2_ = bondAtoms2[i];
        bondList[i].bondType_ = bondMDTypeIndices[i];
    }

    bondFFList = std::vector<CFunctPoint>(state->mdFFBondList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFBondList_.size(); i++)
    {
        CMDFFBond* forceField = state->mdFFBondList_.get(i);
        std::vector<CFunctPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01, rCutoff, numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            bondFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate angle force field vectors
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    mtStateTools.getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices, bondsAcrossPBC);

    angleList = std::vector<CDevAngle>(angleAtoms1.size());
    for(size_t i=0; i<angleAtoms1.size(); i++)
    {
        if(i >= angleAtoms2.size()) continue;
        if(i >= angleAtoms3.size()) continue;
        if(i >= angleMDTypeIndices.size()) continue;

        angleList[i].atomIndex1_ = angleAtoms1[i];
        angleList[i].atomIndex2_ = angleAtoms2[i];
        angleList[i].atomIndex3_ = angleAtoms3[i];
        angleList[i].angleType_ = angleMDTypeIndices[i];
    }

    angleFFList = std::vector<CFunctPoint>(state->mdFFAngleList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFAngleList_.size(); i++)
    {
        CMDFFAngle* forceField = state->mdFFAngleList_.get(i);
        std::vector<CFunctPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            angleFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }

    // Generate dihedral force field vectors
    std::vector<int> dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices;
    mtStateTools.getAllMDDihedralsInSystem(dihedralAtoms1, dihedralAtoms2, dihedralAtoms3, dihedralAtoms4, dihedralMDTypeIndices, bondsAcrossPBC);

    dihedralList = std::vector<CDevDihedral>(dihedralAtoms1.size());
    for(size_t i=0; i<dihedralAtoms1.size(); i++)
    {
        if(i >= dihedralAtoms2.size()) continue;
        if(i >= dihedralAtoms3.size()) continue;
        if(i >= dihedralAtoms4.size()) continue;
        if(i >= dihedralMDTypeIndices.size()) continue;

        dihedralList[i].atomIndex1_ = dihedralAtoms1[i];
        dihedralList[i].atomIndex2_ = dihedralAtoms2[i];
        dihedralList[i].atomIndex3_ = dihedralAtoms3[i];
        dihedralList[i].atomIndex4_ = dihedralAtoms4[i];
        dihedralList[i].dihedralType_ = dihedralMDTypeIndices[i];
    }

    dihedralFFList = std::vector<CFunctPoint>(state->mdFFDihList_.size() * numPointsInForceProfiles);
    for(int i=0; i<state->mdFFDihList_.size(); i++)
    {
        CMDFFDih* forceField = state->mdFFDihList_.get(i);
        std::vector<CFunctPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.0, float(2.0*M_PI), numPointsInForceProfiles));
        for(size_t j=0; j<plot.size(); j++)
        {
            dihedralFFList[toIndexBonded(i, j, numPointsInForceProfiles)] = plot[j];
        }
    }
}

class CFunctorCalcForce
{
public:
    CFunctorCalcForce(int Natoms, int Nbonds, int Natomtypes, int dim, float Lx, float Ly, float Lz, float cutF);

public:
    void setForceFieldMatrices(const CDevForceFieldMatrices& ffMatrices);
    void operator()(CDevAtom& atom);

private:
    C3DVector calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const;
    C3DVector calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i);
    C3DVector calcForceBond(const C3DVector& r_k, const C3DVector& r_i, const int& bondType);
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
    CDevAtom* atomList_;
    CFunctPoint* nonBondFFMatrix_;
    size_t* nonBondFFMatrixFFCount_;
    CDevBond* bondList_;
    CFunctPoint* bondFFList_;
    CDevAngle* angleList_;
    CFunctPoint* angleFFList_;
    CDevDihedral* dihedralList_;
    CFunctPoint* dihedralFFList_;
    CLastError* lastErrorList_;
};

CFunctorCalcForce::CFunctorCalcForce(int Natoms, int Nbonds, int Natomtypes, int dim, float Lx, float Ly, float Lz, float cutF)
{
    Natomtypes_ = Natomtypes;
    Natoms_ = Natoms;
    Nbonds_ = Nbonds;
    dim_ = dim;
    Lx_ = Lx;
    Ly_ = Ly;
    Lz_ = Lz;
    cutF_ = cutF;
}

void CFunctorCalcForce::setForceFieldMatrices(const CDevForceFieldMatrices& ffMatrices)
{
    atomList_ = raw_pointer_cast(&ffMatrices.atomList_[0]);
    nonBondFFMatrix_ = raw_pointer_cast(&ffMatrices.nonBondFFMatrix_[0]);
    nonBondFFMatrixFFCount_ = raw_pointer_cast(&ffMatrices.nonBondFFMatrixFFCount_[0]);
    bondList_ = raw_pointer_cast(&ffMatrices.bondList_[0]);
    bondFFList_ = raw_pointer_cast(&ffMatrices.bondFFList_[0]);
    angleList_ = raw_pointer_cast(&ffMatrices.angleList_[0]);
    angleFFList_ = raw_pointer_cast(&ffMatrices.angleFFList_[0]);
    dihedralList_ = raw_pointer_cast(&ffMatrices.dihedralList_[0]);
    dihedralFFList_ = raw_pointer_cast(&ffMatrices.dihedralFFList_[0]);
    lastErrorList_ = raw_pointer_cast(&ffMatrices.lastErrorList_[0]);
}

void CFunctorCalcForce::operator()(CDevAtom& atom)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx_, 0.0,   0.0 );
    C3DVector PBCy = C3DVector( 0.0,  Ly_,  0.0 );
    C3DVector PBCz = C3DVector( 0.0,  0.0,  Lz_ );

    // Clear forces from primary image (pi = primary image)
    atom.Fpi_ = C3DVector(0.0, 0.0, 0.0);

    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    // :TODO: Later this will be a loop over the neighbor list only!!!
    // :TODO: Here, Coulomb is combined into short range. Should make sure that only Coulomb go past PBC!!!
    int k = atom.index_;
    lastErrorList_[k].reset();
    for(int i=0; i<Natoms_; i++)
    {
        C3DVector r_k = atomList_[k].r_;
        C3DVector r_i = atomList_[i].r_;

        atom.Fpi_+= calcForceNonBondedOn_r_k(r_k, r_i, k, i);

        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCx, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCx, k, i);

        if(dim_ < 2) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCy, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCy, k, i);

        if(dim_ < 3) continue;
        F+= calcForceNonBondedOn_r_k(r_k, r_i + PBCz, k, i);
        F+= calcForceNonBondedOn_r_k(r_k, r_i - PBCz, k, i);
    }
    F+= atom.Fpi_;

    // Add forces from harmonic bonds on particle k
    for(int j=0; j<Nbonds_; j++)
    {
        int iBondTo = -1;
        if((int)bondList_[j].atomIndex1_ == k)
        {
            iBondTo = (int)bondList_[j].atomIndex2_;
        }
        if((int)bondList_[j].atomIndex2_ == k)
        {
            iBondTo = (int)bondList_[j].atomIndex1_;
        }
        if(bondList_[j].bondType_ != -1)
        {
            C3DVector r_k = atomList_[k].r_;
            C3DVector r_i = atomList_[iBondTo].r_;

            atom.Fpi_+= calcForceBond(r_k, r_i, bondList_[j].bondType_);
        }
    }
    F+= atom.Fpi_;

    if(fabs(F.x_) > cutF_) { F.x_ = ((F.x_ >= 0.0) ? 1.0 : -1.0) * cutF_; lastErrorList_[k].lastWarningCode_ = CLastError::warnForcesWereCut; }
    if(fabs(F.y_) > cutF_) { F.y_ = ((F.y_ >= 0.0) ? 1.0 : -1.0) * cutF_; lastErrorList_[k].lastWarningCode_ = CLastError::warnForcesWereCut; }
    if(fabs(F.z_) > cutF_) { F.z_ = ((F.z_ >= 0.0) ? 1.0 : -1.0) * cutF_; lastErrorList_[k].lastWarningCode_ = CLastError::warnForcesWereCut; }

    atom.F_ = F;
}

C3DVector CFunctorCalcForce::calcNonBondForceCoeffs12(const C3DVector& r1, const C3DVector& r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;

    // Return coeffs that correspond to calculating force on r2 (i.e., dr/dx_2, etc.)
    return (r12 * RInv);
}

C3DVector CFunctorCalcForce::calcForceNonBondedOn_r_k(const C3DVector& r_k, const C3DVector& r_i, const int& k, const int& i)
{
    // We know that
    // F_k=-grad_1 U = (-dU/dr) * (dr/dx_k, dr/dy_k, dr/dz_k)
    // We have (-dU/dr) stored in table form from nonBondFFMatrix_,
    // we just need to retrieve the one stored for (k, i) and interpolate
    // it for r = | r_k - r_i |. We use linear interpolation.
    float fSum = 0.0f;
    size_t ffCount = nonBondFFMatrixFFCount_[toIndexNonBond(atomList_[k].typeIndex_, atomList_[i].typeIndex_, Natomtypes_)];
    for(int ffIndex=0; ffIndex<(int)ffCount; ffIndex++)
    {
        float interpVal = nonBondFFMatrix_[toIndexNonBond(atomList_[k].typeIndex_, atomList_[i].typeIndex_, ffIndex, j, Natomtypes_, Natomtypes_, MAX_FF_PER_ATOMIC_SET)];
        fSum+= ...;
    }

    // Now that we have (-dU/dr) at r, we find the analytic calculation
    // of dr/dx_k, dr/dy_k and dr/dz_k, where r = sqrt((r_x_k - r_x_i)^2
    // + (r_y_k - r_y_i)^2 + (r_z_k - r_z_i)^2) and then calulcate the
    // forces on r_k.
    C3DVector c = calcNonBondForceCoeffs12(r_i, r_k);
    return C3DVector();
}

C3DVector CFunctorCalcForce::calcForceBond(const C3DVector& r_k, const C3DVector& r_i, const int& bondType)
{
    return C3DVector();
}

CVelVerlet::CVelVerlet(CMolTwisterState* state, FILE* stdOut)
{
    P = 1.0 / Conv_press;
    p_eps = 1.0;
    eps = 0.0;
    W = 1.0;
    m_dLastFMax = 0.0;
    Fcut = 1000.0;
    m_bCutF = false;
    tau = 1.0;

    const float rCutoff = 10.0f;
    devVelVerlet_ = new CDevForceFieldMatrices(state, stdOut, rCutoff);
}

CVelVerlet::~CVelVerlet()
{
    if(devVelVerlet_) delete devVelVerlet_;
}

void CVelVerlet::Propagator(int N, int dim, double dt, double Lmax, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF, std::vector<C3DVector> &aFpi, bool bNPT)
{
    if(!bNPT)
    {
        double dt_2 = dt / 2.0;
        for(int k=0; k<N; k++)
        {
            double m_k = aParticles[k].m;
            aParticles[k].p+= aF[k]*dt_2;
            aParticles[k].x+= aParticles[k].p*(dt / m_k);
            aF[k] = CalcParticleForce(k, N, dim, Lmax, Lmax, Lmax, aParticles, aFpi[k]);
            aParticles[k].p+= aF[k]*dt_2;
        }
    }
    
    else
    {
        // Fpi = Force periodic image!!!
        // :TODO: run thrust::transform() on each particle (instead of below for loop), where its functor calculates the forces for each particle.
        // The result vector from transform should be aF[k] for all k. To do this, we need to upload some information to the GPU via the transform
        // call. Hence, this information must be included with all particles in the list (the functor will only see a single particle). This therefore
        // includes the following items:
        // * k
        // * N
        // * dim
        // * Lm, Lm, Lm
        // * Fcut
        // * Array of force profiles
        // * Matrix of non-bonded force assignments
        // * Vector of bonded force assignments
        // * Vector of angular force assignments
        // * Vector of dihedral force assignments
        // NOTE! To save space, try to make a common struct, preloaded to the GPU, which thus should be added as members to the functor passed to thrust::transform(). The paricles will each contain
        // * Array of particles (m, r, p, array of neighbour indices)
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, aFpi));  // Step 3.1
        Prop_p(N, dt, aParticles, aF);                      // Step 3.2
        Prop_r(N, dt, aParticles, aF);                      // Step 3.3
        eps+= (dt * p_eps / W);                             // Step 3.4
        double Lm = pow(GetV(Lmax, true), 1.0/3.0);
        for(int k=0; k<N; k++)
            aF[k] = CalcParticleForce(k, N, dim, Lm, Lm, Lm, aParticles, aFpi[k]);
        Prop_p(N, dt, aParticles, aF);                      // Step 3.5
        p_eps+= ((dt / 2.0) * G_eps(N, aParticles, aFpi));  // Step 3.6
    }
}

void CVelVerlet::Prop_p(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<C3DVector>& aF)
{
    double u_eps = p_eps / W;
    double alpha = (1.0 + 1.0/double(N));
    double parm = alpha*u_eps*dt / 4.0;
    double coeff = CMt::Exp(-parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm); // sinh(parm) / parm;
    // :TODO: This can be parallelized using thrust::transform(), where the same array as for Propegator() is used as input (thus, avoiding unecessary uploads)
    for(int k=0; k<N; k++)
    {
        aParticles[k].p = aParticles[k].p*coeff2 + aF[k]*((dt/2.0)*coeff*coeff3);
    }
}

void CVelVerlet::Prop_r(int N, double dt, std::vector<CParticle3D>& aParticles, const std::vector<C3DVector>&)
{
    double u_eps = p_eps / W;
    double parm = u_eps*dt / 2.0;
    double coeff = CMt::Exp(parm);
    double coeff2 = coeff*coeff;
    double coeff3 = CMt::SinhXoverX(parm) * coeff;
    // :TODO: This can be parallelized using thrust::transform(), where the same array as for Propegator() is used as input (thus, avoiding unecessary uploads)
    for(int k=0; k<N; k++)
    {
        C3DVector u_k = aParticles[k].p * (1.0 / aParticles[k].m);
        aParticles[k].x = aParticles[k].x*coeff2 + u_k*(dt*coeff3);
    }
}

C3DVector CVelVerlet::CalcParticleForce(int k, int N, int dim, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles, C3DVector& Fpi)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx, 0.0, 0.0);
    C3DVector PBCy = C3DVector(0.0,  Ly, 0.0);
    C3DVector PBCz = C3DVector(0.0, 0.0,  Lz);
    
    // Add external forces on particle (Fpi is force from primary image only)
    Fpi = C3DVector(0.0, 0.0, 0.0);
    Fpi+= aFExternal[k].CalcForce(aParticles[k].x);
    
    // Add non-bonded forces to particle, as well as
    // non-bonded forces from first PBC images
    for(int i=0; i<N; i++)
    {
        C3DVector r_k = aParticles[k].x;
        C3DVector r_i = aParticles[i].x;
        
        Fpi+= aFNonBonded[i][k].CalcForce(r_k, r_i);
        
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCx);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCx);
        
        if(dim < 2) continue;
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCy);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCy);
        
        if(dim < 3) continue;
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i + PBCz);
        F+= aFNonBonded[i][k].CalcForce(r_k, r_i - PBCz);
    }
    F+= Fpi;
    
    // Add forces from harmonic bonds on particle k
    for(int j=0; j<(int)aFHarmBond.size(); j++)
    {
        int iBondTo = -1;
        int iPotential = -1;
        if(aFHarmBond[j].m_i == k)
        {
            iBondTo = aFHarmBond[j].m_j;
            iPotential = j;
        }
        if(aFHarmBond[j].m_j == k)
        {
            iBondTo = aFHarmBond[j].m_i;
            iPotential = j;
        }
        if(iPotential != -1)
        {
            C3DVector r_k = aParticles[k].x;
            C3DVector r_i = aParticles[iBondTo].x;
            
            Fpi+= aFHarmBond[iPotential].CalcForce(r_k, r_i, Lx, Ly, Lz);
        }
    }
    F+= Fpi;
    
    if(fabs(F.x_) > Fcut) { F.x_ = ((F.x_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    if(fabs(F.y_) > Fcut) { F.y_ = ((F.y_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    if(fabs(F.z_) > Fcut) { F.z_ = ((F.z_ >= 0.0) ? 1.0 : -1.0) * Fcut; m_bCutF = true; }
    StoreMaxF(F);
    PrintDebugInfoAtCutForces(k, N, Lx, Ly, Lz, aParticles);
    
    return F;
}

double CVelVerlet::G_eps(int N, const std::vector<CParticle3D>& aParticles, const std::vector<C3DVector> &aF)
{
    double V = V0 * exp(3.0 * eps);
    
    double sum_p = 0.0;
    double sum_f = 0.0;
    for(int k=0; k<N; k++)
    {
        double p2 = aParticles[k].p * aParticles[k].p;
        sum_p+= (p2 / aParticles[k].m);
        sum_f+= (aF[k] * aParticles[k].x);
    }

    return (1.0 + 1.0 / double(N))*sum_p + sum_f - 3.0*P*V;
}

void CVelVerlet::SetRandMom(double tau)
{
    double a = 2.0 / double(RAND_MAX);
    
    p_eps = (a*double(rand()) - 1.0) * (W / tau);
    if(p_eps == 0.0) p_eps = (W / tau);
}

double CVelVerlet::GetV(double Lmax, bool bNPT)
{
    if(!bNPT) return Lmax*Lmax*Lmax;
    
    return V0 * exp(3.0 * eps);
}

void CVelVerlet::StoreMaxF(const C3DVector& F)
{
    double fAbs = F.norm();
    
    if(fAbs > m_dLastFMax) m_dLastFMax = fAbs;
}

void CVelVerlet::PrintCutMsgAndReset()
{
    if(m_bCutF)
    {
        COut::Printf("\t****************************************\r\n");
        COut::Printf("\t* Warning! Forces were cut!!!          *\r\n");
        COut::Printf("\t****************************************\r\n");
        m_dLastFMax = 0.0;
        m_bCutF = false;
    }
}

void CVelVerlet::PrintDebugInfoAtCutForces(int k, int N, double Lx, double Ly, double Lz, std::vector<CParticle3D>& aParticles)
{
    C3DVector F;
    C3DVector PBCx = C3DVector( Lx, 0.0, 0.0);
    C3DVector PBCy = C3DVector(0.0,  Ly, 0.0);
    C3DVector PBCz = C3DVector(0.0, 0.0,  Lz);

    if(m_bCutF)
    {
        C3DVector r_k = aParticles[k].x;
        COut::Printf("\r\n\r\n------------------ Center -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        
        COut::Printf("\r\n\r\n------------------ Plus PBCx -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCx;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::Printf("\r\n\r\n------------------ Minus PBCx -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCx;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::Printf("\r\n\r\n------------------ Plus PBCy -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCy;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::Printf("\r\n\r\n------------------ Minus PBCy -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCy;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::Printf("\r\n\r\n------------------ Plus PBCz -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x + PBCz;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::Printf("\r\n\r\n------------------ Minus PBCz -----------------------\r\n");
        COut::Printf("Force on particle %i at (%g, %g, %g), Ftot=(%g, %g, %g), |Ftot|=%g\r\n", k, r_k.x_, r_k.y_, r_k.z_, F.x_, F.y_, F.z_, m_dLastFMax);
        for(int i=0; i<N; i++)
        {
            C3DVector r_i = aParticles[i].x - PBCz;
            
            F = aFNonBonded[i][k].CalcForce(r_k, r_i);
            COut::Printf("\t%-4sFrom index %i at (%g, %g, %g) with F=(%g, %g, %g)\r\n", (F.norm() > 10.0) ? "***" : "   ", i, r_i.x_, r_i.y_, r_i.z_, F.x_, F.y_, F.z_);
        }
        COut::CloseOutputFile();
        exit(0);
    }
}

