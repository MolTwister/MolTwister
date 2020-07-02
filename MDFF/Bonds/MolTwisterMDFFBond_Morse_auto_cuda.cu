#include "MolTwisterMDFFBond_Morse.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFBond_Morse::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFBond::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << D_;
        io << alpha_;
        io << r0_;
    }
    else
    {
        io >> D_;
        io >> alpha_;
        io >> r0_;
    }
}

size_t CMDFFBond_Morse::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;

    D_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    alpha_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    r0_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFBond_Morse::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "D:%.6f alpha:%.6f r0:%.6f", J2cal(D_, convertToCal), alpha_, r0_);
    
    return string;
}

std::string CMDFFBond_Morse::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$D$=%g, $\\alpha$=%g, $r_0$=%g", J2cal(D_, convertToCal), alpha_, r0_);
    
    return string;
}

std::string CMDFFBond_Morse::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "morse %.6f %.6f %.6f     # %s %s", J2cal(D_, convertToCal), alpha_, r0_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    
    return line;
}

std::string CMDFFBond_Morse::getHoomdBlueDef(int, int) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "morse.pair_coeff.set('%s', '%s', D0=%.6f, alpha=%.6f, r0=%.6f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
            D_ / lju.energyUnit(), alpha_ * lju.distanceUnit(), r0_ / lju.distanceUnit());
    
    return line;
}

double CMDFFBond_Morse::calcPotential(C3DVector r1, C3DVector r2) const
{
    C3DVector r12 = r2 - r1;
    double R = r12.norm();
    double deltaR = R - r0_;
    double P = (1.0 - exp(-alpha_*deltaR));
    
    return D_ * P*P;
}

void CMDFFBond_Morse::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector r12 = r2 - r1;
    double R = r12.norm();
    double RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double deltaR = R - r0_;
    double Exp = exp(-alpha_*deltaR);
    double RInv_dUdr_neg = -2.0*D_*alpha_*RInv * Exp*(1.0 - Exp);
    C3DVector F = r12*RInv_dUdr_neg;
    
    f2 = F;
    f1 = F*(-1.0);
}

END_CUDA_COMPATIBLE()
