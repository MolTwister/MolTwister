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
#include <vector>
#include <string.h>
#include <string>
#include "Utilities/3DVector.h"
#include "Utilities/CUDAGeneralizations.h"
#include "Utilities/Serializer.h"

BEGIN_CUDA_COMPATIBLE()

class CDefaultAtomicProperties
{
public:
    class CAtomicProperty
    {
    public:
        CAtomicProperty() { }
        CAtomicProperty(std::string ID, C3DVector color, double sigma=1.0, double RCov=0.8) { ID_ = ID; color_ = color; sigma_ = sigma; RCov_ = RCov; }
        
    public:
        C3DVector color_;
        std::string ID_;
        double sigma_;
        double RCov_;
    };
    
public:
    CDefaultAtomicProperties();
    ~CDefaultAtomicProperties();
    
public:
    void serialize(CSerializer& io, bool saveToStream);
    void getCPKColor(std::string ID, double& r, double& g, double& b) const;
    double getWDWRadius(std::string ID) const;
    double getCovalentRadius(std::string ID) const;
    int identifyAtom(std::string ID) const;
    std::string getRecognizedAtomType(std::string ID) const;

private:
    void initAtomicPropertiesArray();
    
private:
    std::vector<CAtomicProperty> atomicProperties_;
};

END_CUDA_COMPATIBLE()
