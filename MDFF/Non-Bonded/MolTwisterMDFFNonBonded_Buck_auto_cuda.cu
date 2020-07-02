#include <math.h>
#include "MolTwisterMDFFNonBonded_Buck.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFNonBonded_Buck::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFNonBonded::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << A_;
        io << rho_;
        io << C_;
    }
    else
    {
        io >> A_;
        io >> rho_;
        io >> C_;
    }
}

std::pair<bool, size_t> CMDFFNonBonded_Buck::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;
    
    A_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    rho_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    C_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    
    return std::pair<bool, size_t>(true, arg);
}

std::string CMDFFNonBonded_Buck::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "A:%.6f rho:%.6f C:%.6f", J2cal(A_, convertToCal), rho_, J2cal(C_, convertToCal));
    
    return string;
}

std::string CMDFFNonBonded_Buck::getArgumentsLaTeX(bool bConvertToCal) const
{
    char string[100];
    
    sprintf(string, "$A$=%g, $\\rho$=%g, $C$=%g", J2cal(A_, bConvertToCal), rho_, J2cal(C_, bConvertToCal));
    
    return string;
}

std::string CMDFFNonBonded_Buck::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "buck/coul/long %.6f %.6f %.6f     # %s %s", J2cal(A_, convertToCal), rho_, J2cal(C_, convertToCal), atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    
    return line;
}

std::string CMDFFNonBonded_Buck::getHoomdBlueDef(int) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "buck.pair_coeff.set('%s', '%s', A=%.6f, rho=%.6f, C=%.6f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
            A_ / lju.energyUnit(), rho_ / lju.distanceUnit(), C_ / lju.buckinghamC());
    
    return line;
}

std::string CMDFFNonBonded_Buck::molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string) const
{
    if(mixWith.getFFType() != getFFType())
    {
        printf("Warning! cannot mix types %s and %s!", mixWith.getFFType().data(), getFFType().data());
        return "";
    }
    
    CMDFFNonBonded_Buck* mixWithBuck = (CMDFFNonBonded_Buck*)&mixWith;
    
    double dA = sqrt(A_ * mixWithBuck->A_);
    double dRho = (rho_ + mixWithBuck->rho_) / 2.0;
    double dC = sqrt(C_ * mixWithBuck->C_);

    char ret[1024];
    sprintf(ret, "%s %s %s %.6f %.6f %.6f", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), getFFType().data(), dA, dRho, dC);

    return ret;
}

double CMDFFNonBonded_Buck::calcPotential(C3DVector r1, C3DVector r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv4 = RInv2*RInv2;
    double      RInv6 = RInv4*RInv2;
    double      rhoInv = (rho_ == 0.0) ? 1E-10 : 1.0 / rho_;
    
    return A_*exp(-R*rhoInv) - C_*RInv6;
}

void CMDFFNonBonded_Buck::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv4 = RInv2*RInv2;
    double      RInv6 = RInv4*RInv2;
    double      rhoInv = (rho_ == 0.0) ? 1E-10 : 1.0 / rho_;
    double      RInv_dUdr_neg = RInv2 * (A_*R*rhoInv*exp(-R*rhoInv) - 6.0*C_*RInv6);
    C3DVector   F = r12*RInv_dUdr_neg;
    
    f2 = F;
    f1 = F*(-1.0);
}

END_CUDA_COMPATIBLE()
