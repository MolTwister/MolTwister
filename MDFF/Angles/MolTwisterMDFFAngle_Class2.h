#pragma once
#include "MolTwisterMDFFAngle.h"

BEGIN_CUDA_COMPATIBLE()

////////////////////////////////////////////////////////////////////////////////////////////////
// Note
////////////////////////////////////////////////////////////////////////////////////////////////
// The COMPASS Class2 forcefield is defined as (using LaTeX notation)
// E = E_a + E_{bb} + E_{ba},
// where
// E_a = K_2 (\theta - \theta_0)^2 + K_3 (\theta - \theta_0)^3 + K_4 (\theta - \theta_0)^4
// E_{bb} = M (r_{ij} - r_1)(r_{jk} - r_2)
// E_{ba} = N_1 (r_{ij} - r_1)(\theta - \theta_0) + N_2 (r_{jk} - r_2)(\theta - \theta_0)
////////////////////////////////////////////////////////////////////////////////////////////////

class CMDFFAngle_Class2 : public CMDFFAngle
{
public:
    class CEaCoeffs
    {
    public:
        CEaCoeffs() { theta0_ = K2_ = K3_ = K4_ = 0.0; }
        
    public:
        double theta0_;
        double K2_;
        double K3_;
        double K4_;
    };
    
    class CEbbCoeffs
    {
    public:
        CEbbCoeffs() { M_ = r1_ = r2_ = 0.0; }
        
    public:
        double M_;
        double r1_;
        double r2_;
    };

    class CEbaCoeffs
    {
    public:
        CEbaCoeffs() { N1_ = N2_ = r1_ = r2_ = 0.0; }
        
    public:
        double N1_;
        double N2_;
        double r1_;
        double r2_;
    };
    
public:
    CMDFFAngle_Class2() : CMDFFAngle() { angleABCEqualToCBA_ = false; toRad_ = 0.017453292519943295; }
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::string getFFType() const { return "Class2"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 3; }
    virtual std::string getHoomdBlueDef(int index, int angleID) const;
    virtual int getNumHoomdBlueDef() const { return 3; }
    virtual std::shared_ptr<CMDFFAngle> createCopy() const { auto ret = std::shared_ptr<CMDFFAngle>(new CMDFFAngle_Class2); *(CMDFFAngle_Class2*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2, C3DVector r3) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector& f1, C3DVector& f2, C3DVector& f3) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Ea.theta0> <Ea.K2> <Ea.K3> <Ea.K4> <Ebb.M> <Ebb.r1> <Ebb.r2> <Eba.N1> <Eba.N2> <Eba.r1> <Eba.r2>" }; }

    CEaCoeffs getEaCoeffs() const { return Ea_; }
    CEbbCoeffs getEbbCoeffs() const { return Ebb_; }
    CEbaCoeffs getEbaCoeffs() const { return Eba_; }
    
private:
    double calcPotentialHarm(C3DVector r1, C3DVector r2, C3DVector r3, double k, double theta0, int degree) const;
    void calcForcesHarm(C3DVector r1, C3DVector r2, C3DVector r3, double k, double theta0, int degree, C3DVector& f1, C3DVector& f2, C3DVector& f3) const;
    double calcPotentialCartesian(C3DVector r1, C3DVector r2, C3DVector r3, double M, double ra, double rb) const;
    void calcForcesCartesian(C3DVector r1, C3DVector r2, C3DVector r3, double M, double ra, double rb, C3DVector& f1, C3DVector& f2, C3DVector& f3) const;
    double calcPotentialMix(C3DVector r1, C3DVector r2, C3DVector r3, double N1, double N2, double ra, double rb, double theta0) const;
    void calcForcesMix(C3DVector r1, C3DVector r2, C3DVector r3, double N1, double N2, double ra, double rb, double theta0, C3DVector& f1, C3DVector& f2, C3DVector& f3) const;
    
private:
    CEaCoeffs Ea_;
    CEbbCoeffs Ebb_;
    CEbaCoeffs Eba_;
    
private:
    double toRad_;
};

END_CUDA_COMPATIBLE()
