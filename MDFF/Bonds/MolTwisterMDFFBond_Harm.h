#pragma once
#include "MolTwisterMDFFBond.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// Note
////////////////////////////////////////////////////////////////////////////////////////////////
// The value given to MolTwister is k, where U=k * (r - r0)^2. For example, to
// LAMMPS scripts we give k and in Hoomd-Blue we must give 2*k.
////////////////////////////////////////////////////////////////////////////////////////////////

class CMDFFBond_Harm : public CMDFFBond
{
public:
    CMDFFBond_Harm() : CMDFFBond() { k_ = 0.0; r0_ = 0.0; }
    
public:
    virtual std::string getFFType() const { return "Harm"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index, int bondID) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::shared_ptr<CMDFFBond> createCopy() const { auto ret = std::shared_ptr<CMDFFBond>(new CMDFFBond_Harm); *(CMDFFBond_Harm*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Harmonic k> <Harmonic r_0>" }; }

    double getK() const { return k_; }
    double getR0() const { return r0_; }
    
private:
    double k_;
    double r0_;
};
