#pragma once
#include "MolTwisterMDFFDih.h"

class CMDFFDih_Harm : public CMDFFDih
{
public:
    CMDFFDih_Harm() : CMDFFDih() { K_ = 0.0; D_ = 1; N_ = 0; }
    
public:
    virtual std::string getFFType() const { return "Harm"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index, int dihedralID) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::shared_ptr<CMDFFDih> createCopy() const { auto ret = std::shared_ptr<CMDFFDih>(new CMDFFDih_Harm); *(CMDFFDih_Harm*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool bConvertToCal) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Harmonic K> <Harmonic D> <Harmonic N>" }; }

    double getK() const { return K_; }
    int getD() const { return D_; }
    int getN() const { return N_; }
    
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4, C3DVector& f1, C3DVector& f2, C3DVector& f3, C3DVector& f4) const;
    
private:
    double K_;
    int D_;
    int N_;
};
