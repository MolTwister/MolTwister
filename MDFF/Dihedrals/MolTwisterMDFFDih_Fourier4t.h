#pragma once
#include "MolTwisterMDFFDih.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFDih_Fourier4t : public CMDFFDih
{
public:
    CMDFFDih_Fourier4t() : CMDFFDih() { V_[0] = V_[1] = V_[2] = V_[3] = 0.0; }
    
public:
    virtual std::string getFFType() const { return "Fourier4t"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index, int dihedralID) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::shared_ptr<CMDFFDih> createCopy() const { auto ret = std::shared_ptr<CMDFFDih>(new CMDFFDih_Fourier4t); *(CMDFFDih_Fourier4t*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<V1> <V2> <V3> <V4>" }; }

    double getV(int index) const { return V_[index]; }

    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4, C3DVector& f1, C3DVector& f2, C3DVector& f3, C3DVector& f4) const;
    
private:
    double V_[4];
};

END_CUDA_COMPATIBLE()
