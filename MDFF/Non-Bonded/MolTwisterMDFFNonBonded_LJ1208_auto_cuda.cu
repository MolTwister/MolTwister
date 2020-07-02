#include "MolTwisterMDFFNonBonded_LJ1208.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFNonBonded_LJ1208::serialize(std::stringstream& io, bool saveToStream)
{
    CMDFFNonBonded::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << sigma_;
        io << epsilon_;
        io << alpha_;
    }
    else
    {
        io >> sigma_;
        io >> epsilon_;
        io >> alpha_;
    }
}

std::pair<bool, size_t> CMDFFNonBonded_LJ1208::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;
    
    epsilon_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    sigma_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    
    if(CASCIIUtility::getArg(arguments, arg++) == "alpha")
    {
        alpha_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    }
    else arg--;
    
    return std::pair<bool, size_t>(true, arg);
}

std::string CMDFFNonBonded_LJ1208::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "epsilon:%.10f sigma:%.10f alpha:%.10f", J2cal(epsilon_, convertToCal), sigma_, alpha_);
    
    return string;
}

std::string CMDFFNonBonded_LJ1208::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$\\epsilon$=%g, $\\sigma$=%g, $\\alpha$=%g", J2cal(epsilon_, convertToCal), sigma_, alpha_);
    
    return string;
}

std::string CMDFFNonBonded_LJ1208::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "lj1208/cut/coul/long %.10f %.10f     # %s %s", J2cal(epsilon_, convertToCal), sigma_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    if(alpha_ != 1.0)
    {
        printf("Warning [LJ1208 (%s,%s)]: LAMMPS does not support the alpha coefficient of Lennard-Jones potentials (see Hoomd-Blue doc.)!", atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    }
    
    return line;
}

std::string CMDFFNonBonded_LJ1208::getHoomdBlueDef(int) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "lj1208.pair_coeff.set('%s', '%s', epsilon=%.10f, sigma=%.10f, alpha=%.10f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
            epsilon_ / lju.energyUnit(), sigma_ / lju.distanceUnit(), alpha_);
    
    return line;
}

std::string CMDFFNonBonded_LJ1208::molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const
{
    if(mixWith.getFFType() != getFFType())
    {
        printf("Warning! cannot mix types %s and %s!", mixWith.getFFType().data(), getFFType().data());
        return "";
    }
    
    CMDFFNonBonded_LJ1208* mixWithLJ = (CMDFFNonBonded_LJ1208*)&mixWith;
    
    double epsilon = sqrt(epsilon_ * mixWithLJ->epsilon_);
    
    double sigma;
    if(mixingRule == "aritmetic")       sigma = (sigma_ + mixWithLJ->sigma_) / 2.0;
    else if(mixingRule == "geometric")  sigma = sqrt(sigma_ * mixWithLJ->sigma_);
    else
    {
        printf("Error! unknown mixing type %s!", mixingRule.data());
        return "";
    }
    
    char ret[1024];
    sprintf(ret, "%s %s %s %.10f %.10f", atomNamesToBond_[0].data(), mixWithLJ->atomNamesToBond_[1].data(), getFFType().data(), epsilon, sigma);

    return ret;
}

double CMDFFNonBonded_LJ1208::calcPotential(C3DVector r1, C3DVector r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv4 = RInv2*RInv2;
    double      RInv6 = RInv4*RInv2;
    double      RInv8 = RInv4*RInv4;
    double      RInv12 = RInv6*RInv6;
    double      sigma2 = sigma_*sigma_;
    double      sigma4 = sigma2*sigma2;
    double      sigma6 = sigma4*sigma2;
    double      sigma8 = sigma4*sigma4;
    double      sigma12 = sigma6*sigma6;
    
    return 4.0*epsilon_ * (sigma12*RInv12 - alpha_*sigma8*RInv8);
}

void CMDFFNonBonded_LJ1208::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv4 = RInv2*RInv2;
    double      RInv6 = RInv4*RInv2;
    double      RInv8 = RInv4*RInv4;
    double      RInv12 = RInv6*RInv6;
    double      sigma2 = sigma_*sigma_;
    double      sigma4 = sigma2*sigma2;
    double      sigma6 = sigma4*sigma2;
    double      sigma8 = sigma4*sigma4;
    double      sigma12 = sigma6*sigma6;
    double      RInv_dUdr_neg = 4.0*epsilon_*RInv2 * (12.0*sigma12*RInv12 - 8.0*alpha_*sigma8*RInv8);
    C3DVector   F = r12*RInv_dUdr_neg;

    f2 = F;
    f1 = F*(-1.0);
}

END_CUDA_COMPATIBLE()
