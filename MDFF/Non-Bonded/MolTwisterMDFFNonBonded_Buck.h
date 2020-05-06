#pragma once
#include "MolTwisterMDFFNonBonded.h"

class CMDFFNonBonded_Buck : public CMDFFNonBonded
{
public:
    CMDFFNonBonded_Buck() : CMDFFNonBonded() { A_ = 0.0; rho_ = 0.0; C_ = 0.0; }
    
public:
    virtual std::string getFFType() const { return "Buck"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::string molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const;
    virtual std::shared_ptr<CMDFFNonBonded> createCopy() const { auto ret = std::shared_ptr<CMDFFNonBonded>(new CMDFFNonBonded_Buck); *(CMDFFNonBonded_Buck*)(ret.get()) = *this; return ret; }
    virtual std::pair<bool, size_t> onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool bConvertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Buckingham A> <Buckingham rho> <Buckingham C>" }; }

    double getA() const { return A_; }
    double getRho() const { return rho_; }
    double getC() const { return C_; }
    
private:
    double A_;
    double rho_;
    double C_;
};
