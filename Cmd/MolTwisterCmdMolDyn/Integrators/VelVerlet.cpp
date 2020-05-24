#include "VelVerlet.h"
#include "../MDLoop/Printf.h"
#include "Math.h"
#include <float.h>
#include <math.h>
#include <functional>

class CDevAtom
{
public:
    CDevAtom() { typeIndex_ = -1; r_[0] = r_[1] = r_[2] = 0.0f; p_[0] = p_[1] = p_[2] = 0.0f; }

public:
    float m_;
    float r_[3];
    float p_[3];
    int typeIndex_;
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

template<class T> T* raw_pointer_cast(T* ptr) { return ptr; }

size_t toIndex(size_t rowIndex, size_t columnIndex, size_t ffIndex, size_t pointIndex, size_t columnCount, size_t rowCount, size_t maxNumFFPerAtomicSet)
{
    return columnCount*(rowCount*(maxNumFFPerAtomicSet*pointIndex + ffIndex) + rowIndex) + columnIndex;
}

size_t toIndex(size_t rowIndex, size_t columnIndex, size_t columnCount)
{
    return columnCount*rowIndex + columnIndex;
}

void prepareFFMatrices(CMolTwisterState* state, std::vector<CDevAtom>& atomList, std::vector<CFunctPoint>& nonBondFFMatrix, std::vector<size_t>& nonBondFFMatrixFFCount)
{
    const float rCutoff = 10.0f;
    const int maxNumFFPerAtomicSet = 5;
    const int numPointsInForceProfiles = 100;

    // Generate atom list for GPU and set initial positions
    size_t numAtoms = state->atoms_.size();
    atomList = std::vector<CDevAtom>(numAtoms);
    for(size_t i=0; i<numAtoms; i++)
    {
        if(state->atoms_[i]->r_.size() == 0) continue;
        atomList[i].r_[0] = state->atoms_[i]->r_[0].x_;
        atomList[i].r_[1] = state->atoms_[i]->r_[0].y_;
        atomList[i].r_[2] = state->atoms_[i]->r_[0].z_;
        atomList[i].m_ = state->atoms_[i]->m_;
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
                nonBondFFMatrixFFCount[toIndex(r, c, numAtomTypes)] = ffIndexList->size();
                for(size_t i=0; i<ffIndexList->size(); i++)
                {
                    int ffIndex = (*ffIndexList)[i];
                    if(ffIndex >= 0)
                    {
                        CMDFFNonBonded* forceField = state->mdFFNonBondedList_.get(ffIndex);
                        std::vector<CFunctPoint> plot = toFuncPtPlot(forceField->calc1DForceProfile(0.01, rCutoff, numPointsInForceProfiles));
                        for(size_t j=0; j<plot.size(); j++)
                        {
                            nonBondFFMatrix[toIndex(r, c, ffIndex, j, numAtomTypes, numAtomTypes, maxNumFFPerAtomicSet)] = plot[j];
                        }
                    }
                }
            }
        }
    }

    // Generate bond force field vector
}

CVelVerlet::CVelVerlet()
{
    P = 1.0 / Conv_press;
    p_eps = 1.0;
    eps = 0.0;
    W = 1.0;
    m_dLastFMax = 0.0;
    Fcut = 1000.0;
    m_bCutF = false;
    tau = 1.0;

    state_ = nullptr;
}

void CVelVerlet::init(CMolTwisterState* state)
{
    state_ = state;
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

void CVelVerlet::Prop_p(int N, double dt, std::vector<CParticle3D>& aParticles, std::vector<C3DVector>& aF)
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

void CVelVerlet::Prop_r(int N, double dt, std::vector<CParticle3D>& aParticles, std::vector<C3DVector>&)
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

double CVelVerlet::G_eps(int N, std::vector<CParticle3D>& aParticles, std::vector<C3DVector> &aF)
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
    double          a = 2.0 / double(RAND_MAX);
    
    p_eps = (a*double(rand()) - 1.0) * (W / tau);
    if(p_eps == 0.0) p_eps = (W / tau);
}

double CVelVerlet::GetV(double Lmax, bool bNPT)
{
    if(!bNPT) return Lmax*Lmax*Lmax;
    
    return V0 * exp(3.0 * eps);
}

void CVelVerlet::StoreMaxF(C3DVector& F)
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

