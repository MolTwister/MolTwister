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

#pragma once
#include "MolTwisterMDFFBond.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFBond_LJC : public CMDFFBond
{
public:
    CMDFFBond_LJC() : CMDFFBond() { epsilon_ = 0.0; sigma_ = 0.0; scale_ = 0.0; }
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::string getFFType() const { return "LJC"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index, int bondID) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::shared_ptr<CMDFFBond> createCopy() const { auto ret = std::shared_ptr<CMDFFBond>(new CMDFFBond_LJC); *(CMDFFBond_LJC*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Lennard-Jones epsilon> <Lennard-Jones sigma> <scale>" }; }

    double getEpsilon() const { return epsilon_; }
    double getSigma() const { return sigma_; }
    double getScale() const { return scale_; }
    
private:
    double epsilon_;
    double sigma_;
    double scale_;
};

END_CUDA_COMPATIBLE()
