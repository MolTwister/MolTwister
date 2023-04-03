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
#include "MolTwisterMDFFAngle.h"

BEGIN_CUDA_COMPATIBLE()

////////////////////////////////////////////////////////////////////////////////////////////////
// Note
////////////////////////////////////////////////////////////////////////////////////////////////
// The value given to MolTwister is k, where U=k * (theta - theta0)^2. For example, to
// LAMMPS scripts we give k and in Hoomd-Blue we must give 2*k.
////////////////////////////////////////////////////////////////////////////////////////////////

class CMDFFAngle_Harm : public CMDFFAngle
{
public:
    CMDFFAngle_Harm() : CMDFFAngle() { k_ = 0.0; theta0_ = 0.0; toRad_ = 0.017453292519943295; }
    
public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::string getFFType() const { return "Harm"; }
    virtual std::string getLammpsDef(int index, bool convertToCal) const;
    virtual int getNumLammpsDef() const { return 1; }
    virtual std::string getHoomdBlueDef(int index, int angleID) const;
    virtual int getNumHoomdBlueDef() const { return 1; }
    virtual std::shared_ptr<CMDFFAngle> createCopy() const { auto ret = std::shared_ptr<CMDFFAngle>(new CMDFFAngle_Harm); *(CMDFFAngle_Harm*)(ret.get()) = *this; return ret; }
    virtual size_t onParse(std::vector<std::string> arguments);
    virtual std::string getArguments(bool convertToCal) const;
    virtual std::string getArgumentsLaTeX(bool convertToCal) const;
    virtual bool isForceCalcAvailable() const { return true; }
    virtual bool isPotentialCalcAvailable() const { return true; }
    virtual double calcPotential(C3DVector r1, C3DVector r2, C3DVector r3) const;
    virtual void calcForces(C3DVector r1, C3DVector r2, C3DVector r3, C3DVector& f1, C3DVector& f2, C3DVector& f3) const;
    virtual std::vector<std::string> getCmdHelpLines() { return { "<Harmonic k> <Harmonic theta0(deg)>" }; }

    double getK() const { return k_; }
    double getTheta0() const { return theta0_; }
    
private:
    double k_;
    double theta0_;

private:
    double toRad_;
};

END_CUDA_COMPATIBLE()
