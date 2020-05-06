#pragma once
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "Utilities/3DVector.h"

////////////////////////////////////////////////////////////////////////////////////////////////
// Note
////////////////////////////////////////////////////////////////////////////////////////////////
// Force is the force on r2 mediated by r1
////////////////////////////////////////////////////////////////////////////////////////////////

class CMDFFBond
{
public:
    enum EDetCrit { critAllRLessThan_R0=0, critMolRLessThan_R0=1, critOnlyVisibleBonds=2, critOnly14Bonds=3 };
    
public:
    CMDFFBond() { bondDetectionCriteria_ = critAllRLessThan_R0; detCritR0_ = 1.2; }
    virtual ~CMDFFBond() {}
    
public:
    virtual std::string getFFType() const = 0;
    virtual std::string getLammpsDef(int index, bool convertToCal) const = 0;
    virtual int getNumLammpsDef() const = 0;
    virtual std::string getHoomdBlueDef(int index, int bondID) const = 0;
    virtual int getNumHoomdBlueDef() const = 0;
    virtual std::shared_ptr<CMDFFBond> createCopy() const = 0;
    void parse(std::vector<std::string> arguments);
    virtual size_t onParse(std::vector<std::string> arguments) = 0;
    virtual std::string getArguments(bool convertToCal) const = 0;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const = 0;
    std::string getAtomInBond(int index) const { return atomNamesToBond_[index]; }
    void setAtomsToBond(std::string atom1, std::string atom2) { atomNamesToBond_[0] = atom1; atomNamesToBond_[1] = atom2; }
    void setBondDetectionCriteria(EDetCrit detCrit, double R0) { bondDetectionCriteria_ = detCrit; detCritR0_ = R0; }
    EDetCrit getBondDetectionCriteria(double& R0) const { R0 = detCritR0_; return bondDetectionCriteria_; }
    std::string getComments() const { return comments_; }
    virtual bool isForceCalcAvailable() const { return false; }
    virtual bool isPotentialCalcAvailable() const { return false; }
    virtual double calcPotential(C3DVector, C3DVector) const { return 0.0; }
    virtual void calcForces(C3DVector, C3DVector, C3DVector& f1, C3DVector& f2) const { f1 = f2 = C3DVector(0.0, 0.0, 0.0); }
    bool calcForcesNumerically(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2, double error=1E-3, int maxIter=50) const;
    virtual std::vector<std::string> getCmdHelpLines() = 0;

protected:
    double J2cal(double val, bool convertToCal) const { return convertToCal ? (val / 4.184) : val; }
    
protected:
    std::string atomNamesToBond_[2];
    std::string comments_;
    
private:
    EDetCrit bondDetectionCriteria_;
    double detCritR0_;
};
