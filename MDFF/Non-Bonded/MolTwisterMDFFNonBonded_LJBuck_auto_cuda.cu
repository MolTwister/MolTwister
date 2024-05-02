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

#include "MolTwisterMDFFNonBonded_LJBuck.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"
#include <math.h>

BEGIN_CUDA_COMPATIBLE()

void CMDFFNonBonded_LJBuck::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFNonBonded::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << ljBuckType_;
        io << sigma_;
        io << epsilon_;
        io << A_;
        io << rho_;
        io << C_;
    }
    else
    {
        io >> ljBuckType_;
        io >> sigma_;
        io >> epsilon_;
        io >> A_;
        io >> rho_;
        io >> C_;
    }
}

std::pair<bool, size_t> CMDFFNonBonded_LJBuck::onParse(std::vector<std::string> arguments)
{
    size_t arg = 0;

    ljBuckType_ = atoi(CASCIIUtility::getArg(arguments, arg++).data());
    
    if(ljBuckType_ == 0)
    {
        epsilon_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
        sigma_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    }
    else if(ljBuckType_ == 1)
    {
        A_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
        rho_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
        C_ = atof(CASCIIUtility::getArg(arguments, arg++).data());
    }
    else
    {
        printf("Error: Invalid LJBuck type %i!\r\n", ljBuckType_);
        return std::pair<bool, size_t>(false, arg);
    }
    
    return std::pair<bool, size_t>(true, arg);
}

std::string CMDFFNonBonded_LJBuck::getArguments(bool convertToCal) const
{
    char string[100] = "NA";
    
    if(ljBuckType_ == 0)
        sprintf(string, "epsilon:%.10f sigma:%.10f", J2cal(epsilon_, convertToCal), sigma_);
    else if(ljBuckType_ == 1)
        sprintf(string, "A:%.6f rho:%.6f C:%.6f", J2cal(A_, convertToCal), rho_, J2cal(C_, convertToCal));
    
    return string;
}

std::string CMDFFNonBonded_LJBuck::getArgumentsLaTeX(bool convertToCal) const
{
    char string[100] = "NA";
    
    if(ljBuckType_ == 0)
        sprintf(string, "$\\epsilon$=%g, $\\sigma$=%g", J2cal(epsilon_, convertToCal), sigma_);
    else if(ljBuckType_ == 1)
        sprintf(string, "$A$=%g, $\\rho$=%g, $C$=%g", J2cal(A_, convertToCal), rho_, J2cal(C_, convertToCal));
    
    return string;
}

std::string CMDFFNonBonded_LJBuck::getLammpsDef(int, bool convertToCal) const
{
    char line[1024] = "Error unknown LJBuck type ID!";
    
    if(ljBuckType_ == 0)
        sprintf(line, "ljbuck/cut/coul/long 0 %.10f %.10f     # %s %s", J2cal(epsilon_, convertToCal), sigma_, atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    else if(ljBuckType_ == 1)
        sprintf(line, "ljbuck/cut/coul/long 1 %.6f %.6f %.6f     # %s %s", J2cal(A_, convertToCal), rho_, J2cal(C_, convertToCal), atomNamesToBond_[0].data(), atomNamesToBond_[1].data());
    
    return line;
}

std::string CMDFFNonBonded_LJBuck::getHoomdBlueDef(int) const
{
    char line[1024] = "Error unknown LJBuck type ID!";
    CLJUnits lju;
    
    if(ljBuckType_ == 0)
    {
        sprintf(line, "lj.pair_coeff.set('%s', '%s', epsilon=%.10f, sigma=%.10f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
                epsilon_ / lju.energyUnit(), sigma_ / lju.distanceUnit());
    }
    else if(ljBuckType_ == 1)
    {
        sprintf(line, "buck.pair_coeff.set('%s', '%s', A=%.6f, rho=%.6f, C=%.6f)", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(),
                A_ / lju.energyUnit(), rho_ / lju.distanceUnit(), C_ / lju.buckinghamC());
    }
    
    return line;
}

std::string CMDFFNonBonded_LJBuck::molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const
{
    if(mixWith.getFFType() != getFFType())
    {
        printf("Warning! cannot mix types %s and %s!", mixWith.getFFType().data(), getFFType().data());
        return "";
    }

    CMDFFNonBonded_LJBuck* mixWithLJBuck = (CMDFFNonBonded_LJBuck*)&mixWith;
    if(mixWithLJBuck->ljBuckType_ != ljBuckType_)
    {
        printf("Error! cannot mix Lennard-Jones and Buckingham parameters!");
        return "";
    }
    
    char ret[1024];
    if(ljBuckType_ == 0)
    {
        double epsilon = sqrt(epsilon_ * mixWithLJBuck->epsilon_);
        
        double sigma;
        if(mixingRule == "aritmetic")       sigma = (sigma_ + mixWithLJBuck->sigma_) / 2.0;
        else if(mixingRule == "geometric")  sigma = sqrt(sigma_ * mixWithLJBuck->sigma_);
        else
        {
            printf("Error! unknown mixing type %s!", mixingRule.data());
            return "";
        }
        
        sprintf(ret, "%s %s %s %.10f %.10f", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), getFFType().data(), epsilon, sigma);
    }
    else if(ljBuckType_ == 1)
    {
        double A = sqrt(A_ * mixWithLJBuck->A_);
        double rho = (rho_ + mixWithLJBuck->rho_) / 2.0;
        double C = sqrt(C_ * mixWithLJBuck->C_);
        
        sprintf(ret, "%s %s %s %.6f %.6f %.6f", atomNamesToBond_[0].data(), atomNamesToBond_[1].data(), getFFType().data(), A, rho, C);
    }
    else
    {
        printf("Error! unknown LJBuck type %i!", ljBuckType_);
        return "";
    }

    return ret;
}

double CMDFFNonBonded_LJBuck::calcPotential(C3DVector r1, C3DVector r2) const
{
    if(ljBuckType_ == 0)
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
        
        return 4.0*epsilon_ * (sigma12*RInv12 - sigma6*RInv6);
    }

    if(ljBuckType_ == 1)
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
    
    return 0.0;
}

void CMDFFNonBonded_LJBuck::calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const
{
    if(ljBuckType_ == 0)
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
        double      RInv_dUdr_neg = 4.0*epsilon_*RInv2 * (12.0*sigma12*RInv12 - 6.0*sigma6*RInv6);
        C3DVector   F = r12*RInv_dUdr_neg;
        
        f2 = F;
        f1 = F*(-1.0);
    }
    
    if(ljBuckType_ == 1)
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
    
    f1 = f2 = C3DVector(0.0, 0.0, 0.0);
}

END_CUDA_COMPATIBLE()
