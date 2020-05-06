#include "MolTwisterMDFFBond_LJC.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"

size_t CMDFFBond_LJC::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;

    epsilon_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    sigma_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    scale_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFBond_LJC::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "epsilon:%.10f sigma:%.10f scale:%.4f", J2cal(epsilon_, convertToCal), sigma_, scale_);
    
    return string;
}

std::string CMDFFBond_LJC::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$\\epsilon$=%g, $\\sigma$=%g, $\\eta$=%g", J2cal(epsilon_, convertToCal), sigma_, scale_);
    
    return string;
}

std::string CMDFFBond_LJC::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "ljcbond %.6f %.6f %.4f     # %s %s", J2cal(epsilon_, convertToCal), sigma_, scale_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    
    return line;
}

std::string CMDFFBond_LJC::getHoomdBlueDef(int, int bondID) const
{
    char line[1024];
    CLJUnits lju;
    
    sprintf(line, "ljcbond.bond_coeff.set('bondtype%i', epsilon=%.6f, sigma=%.6f, scale=%.4f)", bondID, 2.0*epsilon_ / lju.energyUnit(), sigma_ / lju.distanceUnit(), scale_);
    
    return line;
}
