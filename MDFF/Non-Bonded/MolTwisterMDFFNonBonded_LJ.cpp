//
// Copyright (C) 2021 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include "MolTwisterMDFFNonBonded_LJ.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFNonBonded_LJ::serialize(CSerializer& io, bool saveToStream)
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

std::pair<bool, size_t> CMDFFNonBonded_LJ::onParse(std::vector<std::string> arguments)
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

std::string CMDFFNonBonded_LJ::getArguments(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "epsilon:%.10f sigma:%.10f alpha:%.10f", J2cal(epsilon_, convertToCal), sigma_, alpha_);
    
    return string;
}

std::string CMDFFNonBonded_LJ::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100];
    
    sprintf(string, "$\\epsilon$=%g, $\\sigma$=%g, $\\alpha$=%g", J2cal(epsilon_, convertToCal), sigma_, alpha_);
    
    return string;
}

std::string CMDFFNonBonded_LJ::getLammpsDef(int, bool convertToCal) const
{
    char line[1024];
    
    sprintf(line, "lj/cut/coul/long %.10f %.10f     # %s %s", J2cal(epsilon_, convertToCal), sigma_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    if(alpha_ != 1.0)
    {
        printf("Warning [LJ (%s,%s)]: LAMMPS does not support the alpha coefficient of Lennard-Jones potentials (see Hoomd-Blue doc.)!", atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    }

    return line;
}

std::string CMDFFNonBonded_LJ::getHoomdBlueDef(int) const
{
    char line[1024];
    CLJUnits lju;

    sprintf(line, "lj.pair_coeff.set('%s', '%s', epsilon=%.10f, sigma=%.10f, alpha=%.10f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
            epsilon_ / lju.energyUnit(), sigma_ / lju.distanceUnit(), alpha_);
    
    return line;
}

std::string CMDFFNonBonded_LJ::molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const
{
    if(mixWith.getFFType() != getFFType())
    {
        printf("Warning! cannot mix types %s and %s!", mixWith.getFFType().data(), getFFType().data());
        return "";
    }
    
    CMDFFNonBonded_LJ* mixWithLJ = (CMDFFNonBonded_LJ*)&mixWith;
    
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

double CMDFFNonBonded_LJ::calcPotential(C3DVector r1, C3DVector r2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv6 = RInv2*RInv2*RInv2;
    double      RInv12 = RInv6*RInv6;
    double      sigma2 = sigma_*sigma_;
    double      sigma6 = sigma2*sigma2*sigma2;
    double      sigma12 = sigma6*sigma6;
    
    return 4.0*epsilon_ * (sigma12*RInv12 - alpha_*sigma6*RInv6);
}

void CMDFFNonBonded_LJ::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    C3DVector   r12 = r2 - r1;
    double      R = r12.norm();
    double      RInv = (R == 0.0) ? 1.0 / 1E-10 : 1.0 / R;
    double      RInv2 = RInv*RInv;
    double      RInv6 = RInv2*RInv2*RInv2;
    double      RInv12 = RInv6*RInv6;
    double      sigma2 = sigma_*sigma_;
    double      sigma6 = sigma2*sigma2*sigma2;
    double      sigma12 = sigma6*sigma6;
    double      RInv_dUdr_neg = 4.0*epsilon_*RInv2 * (12.0*sigma12*RInv12 - 6.0*alpha_*sigma6*RInv6);
    C3DVector   F = r12*RInv_dUdr_neg;

    f2 = F;
    f1 = F*(-1.0);
}

END_CUDA_COMPATIBLE()
