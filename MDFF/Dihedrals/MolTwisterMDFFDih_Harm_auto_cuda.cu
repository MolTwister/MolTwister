#include "MolTwisterMDFFDih_Harm.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

size_t CMDFFDih_Harm::onParse(std::vector<std::string> arguments)
{    
    size_t arg = 0;

    K_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    D_ = atoi(CASCIIUtility::getArg(arguments, arg++).data());
    N_ = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFDih_Harm::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "k:%.6f d:%i n:%i", J2cal(K_, convertToCal), D_, N_);
    
    return string;
}

std::string CMDFFDih_Harm::getArgumentsLaTeX(bool bConvertToCal) const
{
    char string[100];
    
    sprintf(string, "$k$=%g, $d$=%i, $n$=%i", J2cal(K_, bConvertToCal), D_, N_);
    
    return string;
}

std::string CMDFFDih_Harm::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "%.6f %i %i     # %s %s %s %s", J2cal(K_, convertToCal), D_, N_,
            atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data(), atomNamesToBond_[3].data());
    
    return line;
}

std::string CMDFFDih_Harm::getHoomdBlueDef(int, int dihedralID) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "dihharm.set_coeff('dihedraltype%i', k=%.6f, d=%i, n=%i)", dihedralID,
            2.0*K_ / lju.energyUnit(), D_, N_);
    
    return line;
}

double CMDFFDih_Harm::calcPotential(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    C3DVector       n1 = r23.cross(r21);
    C3DVector       n2 = r32.cross(r34);
    double          n1n2 = n1.norm()*n2.norm();
    double          cosPhi = (n1n2 == 0.0) ? (n1*n2) / 1E-20 : (n1*n2) / n1n2;
    double          phi = acos(cosPhi);
    
    return K_ * (1.0 + double(D_)*cos(double(N_)*phi));
}

void CMDFFDih_Harm::calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4, C3DVector& f1, C3DVector& f2, C3DVector& f3, C3DVector& f4) const
{
    C3DVector       r21 = r1 - r2;
    C3DVector       r23 = r3 - r2;
    C3DVector       r32 = r23*(-1.0);
    C3DVector       r34 = r4 - r3;
    C3DVector       n1 = r23.cross(r21);
    C3DVector       n2 = r32.cross(r34);
    double          n1n2 = n1.norm()*n2.norm();
    double          cosPhi = (n1n2 == 0.0) ? (n1*n2) / 1E-20 : (n1*n2) / n1n2;
    double          phi = acos(cosPhi);
    double          dU_dPhi = -K_*double(D_*N_) * sin(double(N_)*phi);
    C3DVector       dPhi_dr1 = calcDihedralForceCoeffs14(r4, r3, r2, r1);
    C3DVector       dPhi_dr2 = calcDihedralForceCoeffs23(r4, r3, r2, r1);
    C3DVector       dPhi_dr3 = calcDihedralForceCoeffs23(r1, r2, r3, r4);
    C3DVector       dPhi_dr4 = calcDihedralForceCoeffs14(r1, r2, r3, r4);
    
    f1 = dPhi_dr1*dU_dPhi;
    f2 = dPhi_dr2*dU_dPhi;
    f3 = dPhi_dr3*dU_dPhi;
    f4 = dPhi_dr4*dU_dPhi;
}

END_CUDA_COMPATIBLE()
