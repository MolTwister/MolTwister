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

#pragma once
#include "MolTwisterMDFFNonBonded.h"

BEGIN_CUDA_COMPATIBLE()

class CMDFFNonBonded_LJBuck : public CMDFFNonBonded
{
public:
    CMDFFNonBonded_LJBuck() : CMDFFNonBonded() { sigma_ = 0.0; epsilon_ = 0.0; A_ = 0.0; rho_ = 1.0; C_ = 0.0; ljBuckType_ = 0; }
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::string getFFType() const { return "LJBuck"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::string molTwisterMixCmdArg(const CMDFFNonBonded& mixWith, std::string mixingRule) const;
    virtual std::shared_ptr<CMDFFNonBonded> createCopy() const { auto ret = std::shared_ptr<CMDFFNonBonded>(new CMDFFNonBonded_LJBuck); *(CMDFFNonBonded_LJBuck*)(ret.get()) = *this; return ret; }
    virtual std::pair<bool, size_t> onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector& f1, C3DVector& f2) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "0 <Lennard-Jones epsilon> <Lennard-Jones sigma>", "1 <Buckingham A> <Buckingham rho> <Buckingham C>" }; }

    double getLJBuckType() const { return ljBuckType_; }
    double getSigma() const { return sigma_; }
    double getEpsilon() const { return epsilon_; }
    
private:
    int ljBuckType_; // 0:LJ, 1:Buck
    double sigma_;
    double epsilon_;
    double A_;
    double rho_;
    double C_;
};

END_CUDA_COMPATIBLE()
