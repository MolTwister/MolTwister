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
#include <string>

struct SMolDynConfigStruct
{
    enum Ensemble { ensembleNVT=0, ensembleNPT=1, ensembleNVE=2 };

    Ensemble ensemble_;
    double temperature_;
    double temperatureRelaxationTime_;
    int temperatureNHChainLength_;
    int temperatureRESPASteps_;
    double pressure_;
    double pressureRelaxationTime_;
    int pressureNHChainLength_;
    int pressureRESPASteps_;
    double cutoffRadius_;
    double neighListShell_;
    double cutoffForce_;
    int numberOfTimeSteps_;
    double timeStep_;
    int outputStride_;
    std::string outInfoFile_;
    std::string outXYZFile_;
    std::string outDCDFile_;
    std::string outXTCFile_;
    std::string outPDistrFile_;
    std::string outVDistrFile_;
    bool includeXYZFile_;
    bool includeDCDFile_;
    bool includeXTCFile_;
    float xtcPrecision_;
    bool includePDistrFile_;
    bool includeVDistrFile_;
    double maxPDistrOutput_;
    double maxVDistrOutput_;
    bool verboseOutput_;
    double scale12Interactions_;
    double scale13Interactions_;
    double scale14Interactions_;
    double scaleAbove14BondedInteractions_;
    double gradientDescentLearningRate_;
    int gradientDescentMaxSteps_;
    double gradientDescentAccuracy_;
};
