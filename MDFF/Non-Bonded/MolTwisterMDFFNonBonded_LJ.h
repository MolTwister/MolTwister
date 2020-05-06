#pragma once
#include "MolTwisterMDFFNonBonded.h"

class CMDFFNonBonded_LJ : public CMDFFNonBonded
{
public:
    CMDFFNonBonded_LJ() : CMDFFNonBonded() { sigma_ = 0.0; epsilon_ = 0.0; alpha_ = 1.0; }
    
public:
    virtual std::string getFFType() const { return "LJ"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::string molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const;
    virtual std::shared_ptr<CMDFFNonBonded> createCopy() const { auto ret = std::shared_ptr<CMDFFNonBonded>(new CMDFFNonBonded_LJ); *(CMDFFNonBonded_LJ*)(ret.get()) = *this; return ret; }
    virtual std::pair<bool, size_t> onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Lennard-Jones sigma> <Lennard-Jones epsilon> [alpha <Lennard-Jones alpha>]" }; }

    double getSigma() const { return sigma_; }
    double getEpsilon() const { return epsilon_; }
    
private:
    double sigma_;
    double epsilon_;
    double alpha_;
};
