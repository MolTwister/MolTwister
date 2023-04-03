//
// Copyright (C) 2023 Richard Olsen.
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

#include "MolTwisterMDFFBond_LJC.h"
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/ASCIIUtility.h"

BEGIN_CUDA_COMPATIBLE()

void CMDFFBond_LJC::serialize(CSerializer& io, bool saveToStream)
{
    CMDFFBond::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << epsilon_;
        io << sigma_;
        io << scale_;
    }
    else
    {
        io >> epsilon_;
        io >> sigma_;
        io >> scale_;
    }
}

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

END_CUDA_COMPATIBLE()
