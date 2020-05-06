#include "MolTwisterMDFFDih_Fourier4t.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

size_t CMDFFDih_Fourier4t::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;
    for(int i=0; i<4; i++)
    {
        V_[i] = atof(CASCIIUtility::getArg(arguments, arg++).data());
    }

    return arg;
}

std::string CMDFFDih_Fourier4t::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "V1:%.6f V2:%.6f V3:%.6f V4:%.6f", J2cal(V_[0], convertToCal), J2cal(V_[1], convertToCal), J2cal(V_[2], convertToCal), J2cal(V_[3], convertToCal));
    
    return string;
}

std::string CMDFFDih_Fourier4t::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$V_1$=%g, $V_2$=%g, $V_3$=%g, $V_4$=%g", J2cal(V_[0], convertToCal), J2cal(V_[1], convertToCal), J2cal(V_[2], convertToCal), J2cal(V_[3], convertToCal));
    
    return string;
}

std::string CMDFFDih_Fourier4t::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "%.6f %.6f %.6f %.6f     # %s %s %s %s", J2cal(V_[0], convertToCal), J2cal(V_[1], convertToCal), J2cal(V_[2], convertToCal), J2cal(V_[3], convertToCal),
            atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), atomNamesToBond_[2].data(), atomNamesToBond_[3].data());
    
    return line;
}

std::string CMDFFDih_Fourier4t::getHoomdBlueDef(int, int dihedralID) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "dihfour.set_coeff('dihedraltype%i', k1=%.6f, k2=%.6f, k3=%.6f, k4=%.6f)", dihedralID,
            V_[0] / lju.energyUnit(), V_[1] / lju.energyUnit(), V_[2] / lju.energyUnit(), V_[3] / lju.energyUnit());
    
    return line;
}

double CMDFFDih_Fourier4t::calcPotential(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4) const
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
    double          U1 = V_[0] * (1.0 + cos(phi));
    double          U2 = V_[1] * (1.0 - cos(2.0*phi));
    double          U3 = V_[2] * (1.0 + cos(3.0*phi));
    double          U4 = V_[3] * (1.0 - cos(4.0*phi));
    
    return 0.5 * (U1 + U2 + U3 + U4);
}

void CMDFFDih_Fourier4t::calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector r4, C3DVector& f1, C3DVector& f2, C3DVector& f3, C3DVector& f4) const
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
    double          dU_dPhi1 = -0.5*V_[0]*sin(phi);
    double          dU_dPhi2 = V_[1]*sin(2.0*phi);
    double          dU_dPhi3 = -1.5*V_[2]*sin(3.0*phi);
    double          dU_dPhi4 = 2.0*V_[3]*sin(4.0*phi);
    double          dU_dPhi = dU_dPhi1 + dU_dPhi2 + dU_dPhi3 + dU_dPhi4;
    C3DVector       dPhi_dr1 = calcDihedralForceCoeffs14(r4, r3, r2, r1);
    C3DVector       dPhi_dr2 = calcDihedralForceCoeffs23(r4, r3, r2, r1);
    C3DVector       dPhi_dr3 = calcDihedralForceCoeffs23(r1, r2, r3, r4);
    C3DVector       dPhi_dr4 = calcDihedralForceCoeffs14(r1, r2, r3, r4);
    
    f1 = dPhi_dr1*dU_dPhi;
    f2 = dPhi_dr2*dU_dPhi;
    f3 = dPhi_dr3*dU_dPhi;
    f4 = dPhi_dr4*dU_dPhi;
}
