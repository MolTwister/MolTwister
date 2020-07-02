#include "MolTwisterMDFFBond_Harm.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"

BEGIN_CUDA_COMPATIBLE()

void CMDFFBond_Harm::serialize(std::stringstream& io, bool saveToStream)
{
    CMDFFBond::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << k_;
        io << r0_;
    }
    else
    {
        io >> k_;
        io >> r0_;
    }
}

size_t CMDFFBond_Harm::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;
    
    k_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    r0_ = atof(CASCIIUtility::getArg(arguments, arg++).data());

    return arg;
}

std::string CMDFFBond_Harm::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "k:%.6f r0:%.6f", J2cal(k_, convertToCal), r0_);
    
    return string;
}

std::string CMDFFBond_Harm::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$k$=%g, $r_0$=%g", J2cal(k_, convertToCal), r0_);
    
    return string;
}

std::string CMDFFBond_Harm::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "harmonic %.6f %.6f     # %s %s", J2cal(k_, convertToCal), r0_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    
    return line;
}

std::string CMDFFBond_Harm::getHoomdBlueDef(int, int bondID) const
{
    char line[1024];
    CLJUnits lju;
  
    sprintf(line, "harmbond.bond_coeff.set('bondtype%i', k=%.6f, r0=%.6f)", bondID, 2.0*k_ / lju.harmonicBondK(), r0_ / lju.distanceUnit());
    
    return line;
}

double CMDFFBond_Harm::calcPotential(C3DVector r1, C3DVector r2) const
{
    C3DVector r12 = r2 - r1;
    double R = r12.norm();
    double deltaR = R - r0_;
    
    return k_ * deltaR*deltaR;
}

void CMDFFBond_Harm::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector r12 = r2 - r1;
    double R = r12.norm();
    double RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double deltaR = R - r0_;
    double RInv_dUdr_neg = -2.0*k_ * deltaR*RInv;
    C3DVector F = r12*RInv_dUdr_neg;
    
    f2 = F;
    f1 = F*(-1.0);
}

END_CUDA_COMPATIBLE()
