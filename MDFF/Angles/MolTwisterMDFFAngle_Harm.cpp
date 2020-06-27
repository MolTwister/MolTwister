#include <math.h>
#include "MolTwisterMDFFAngle_Harm.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

size_t CMDFFAngle_Harm::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;

    k_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    theta0_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFAngle_Harm::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "k:%.6f theta0:%.6f", J2cal(k_, convertToCal), theta0_);
    
    return string;
}

std::string CMDFFAngle_Harm::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$k$=%g, $\\theta_0$=%g", J2cal(k_, convertToCal), theta0_);
    
    return string;
}

std::string CMDFFAngle_Harm::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "harmonic %.6f %.6f     # %s %s %s", J2cal(k_, convertToCal), theta0_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data());
    
    return line;
}

std::string CMDFFAngle_Harm::getHoomdBlueDef(int, int angleID) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "harmangle.set_coeff('angletype%i', k=%.6f, t0=%.6f)", angleID, 2.0*k_ / lju.energyUnit(), theta0_ * M_PI/180.0);
    
    return line;
}

double CMDFFAngle_Harm::calcPotential(C3DVector r1, C3DVector r2, C3DVector r3) const
{
    C3DVector   r21 = r1 - r2;
    C3DVector   r23 = r3 - r2;
    double      R21 = r21.norm();
    double      R23 = r23.norm();
    double      R21Inv = (R21 == 0.0) ? 1.0 / 1E-10 : 1.0 / R21;
    double      R23Inv = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      cosTheta = (r21*r23) * (R21Inv*R23Inv);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0_ * toRad_;
    double      deltaTheta = Theta - Theta0;
    
    return k_ * deltaTheta*deltaTheta;
}

void CMDFFAngle_Harm::calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector& f1, C3DVector& f2, C3DVector& f3) const
{
    C3DVector   r21 = r1 - r2;
    C3DVector   r23 = r3 - r2;
    double      R21 = r21.norm();
    double      R23 = r23.norm();
    double      R21Inv = (R21 == 0.0) ? 1.0 / 1E-10 : 1.0 / R21;
    double      R23Inv = (R23 == 0.0) ? 1.0 / 1E-10 : 1.0 / R23;
    double      cosTheta = (r21*r23) * (R21Inv*R23Inv);
    double      Theta = acos(cosTheta);
    double      Theta0 = theta0_ * toRad_;
    double      deltaTheta = Theta - Theta0;
    C3DVector   dTheta_dr3 = calcAngularForceCoeffs13(r1, r2, r3);
    C3DVector   dTheta_dr2 = calcAngularForceCoeffs2(r1, r2, r3);
    C3DVector   dTheta_dr1 = calcAngularForceCoeffs13(r3, r2, r1);
    double      dU_dTheta = -2.0*k_ * deltaTheta;

    f3 = dTheta_dr3*dU_dTheta;
    f2 = dTheta_dr2*dU_dTheta;
    f1 = dTheta_dr1*dU_dTheta;
}

END_CUDA_COMPATIBLE()
