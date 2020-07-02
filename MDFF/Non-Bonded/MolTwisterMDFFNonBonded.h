#pragma once
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <sstream>
#include "Utilities/3DVector.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

////////////////////////////////////////////////////////////////////////////////////////////////
// Note
////////////////////////////////////////////////////////////////////////////////////////////////
// Force is the force on r2 mediated by r1 if we define r_vec = r2_vec - r1_vec
////////////////////////////////////////////////////////////////////////////////////////////////

class CMDFFNonBonded
{
public:
    CMDFFNonBonded() {}
    virtual ~CMDFFNonBonded() {}
    
public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::string getFFType() const = 0;
    virtual std::string getLammpsDef(int index, bool convertToCal) const = 0;
    virtual int getNumLammpsDef() const = 0;
    virtual std::string getHoomdBlueDef(int index) const = 0;
    virtual int getNumHoomdBlueDef() const = 0;
    virtual std::string molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const = 0;
    virtual std::shared_ptr<CMDFFNonBonded> createCopy() const = 0;
    bool parse(std::vector<std::string> arguments);
    virtual std::pair<bool, size_t> onParse(std::vector<std::string> arguments) = 0;
    virtual std::string getArguments(bool convertToCal) const = 0;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const = 0;
    std::string getAtomInBond(int index) const { return atomNamesToBond_[index]; }
    void setAtomsToBond(std::string atom1, std::string atom2) { atomNamesToBond_[0] = atom1; atomNamesToBond_[1] = atom2; }
    std::string getComments() const { return comments_; }
    virtual bool isForceCalcAvailable() const { return false; }
    virtual bool isPotentialCalcAvailable() const { return false; }
    virtual double calcPotential(C3DVector, C3DVector) const { return 0.0; }
    virtual void calcForces(C3DVector, C3DVector, C3DVector& f1, C3DVector& f2) const { f1 = f2 = C3DVector(0.0, 0.0, 0.0); }
    bool calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2, double error=1E-3, int maxIter=50) const;
    virtual std::vector<std::string> getCmdHelpLines() = 0;
    std::vector<std::pair<float, float>> calc1DForceProfile(float rStart, float rEnd, int points) const;
    std::vector<std::pair<float, float>> calc1DPotentialProfile(float rStart, float rEnd, int points) const;

protected:
    double J2cal(double val, bool convertToCal) const { return convertToCal ? (val / 4.184) : val; }
    C3DVector calcNonBondForceCoeffs12(C3DVector r1, C3DVector r2) const;

protected:
    std::string atomNamesToBond_[2];
    std::string comments_;
};

END_CUDA_COMPATIBLE()
