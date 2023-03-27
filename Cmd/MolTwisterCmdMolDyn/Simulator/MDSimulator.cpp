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

#include "MDSimulator.h"
#include "../../../MolTwisterState.h"
#include "../MDLoop/MDLoop.h"
#include "../MDLoop/Printf.h"

BEGIN_CUDA_COMPATIBLE()

void CMDSimulator::run(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent)
{
    CMolTwisterState state;
    state.serialize(stateContent, false);
    run(config, stdOut, &state);
}

void CMDSimulator::run(SMolDynConfigStruct config, FILE* stdOut, void* state)
{
    COut::setStdOut(stdOut);

    CMDLoop mdLoop(config.includeXYZFile_, config.outXYZFile_, config.includeDCDFile_, config.outDCDFile_);
    mdLoop.setPDistrOutput(config.includePDistrFile_, config.maxPDistrOutput_, config.outPDistrFile_);
    mdLoop.setVDistrOutput(config.includeVDistrFile_, config.maxVDistrOutput_, config.outVDistrFile_);
    CSimulationBox simBox((CMolTwisterState*)state, stdOut, config);

    COut::setOutputFile(fopen(config.outInfoFile_.data(), "w"));
    mdLoop.runSimulation(simBox, config.numberOfTimeSteps_, config.outputStride_);
    COut::closeOutputFile();
}

void CMDSimulator::optimize(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent)
{
    CMolTwisterState state;
    state.serialize(stateContent, false);
    optimize(config, stdOut, &state);
}

void CMDSimulator::optimize(SMolDynConfigStruct config, FILE* stdOut, void* state)
{
    COut::setStdOut(stdOut);

    CMDLoop mdLoop(config.includeXYZFile_, config.outXYZFile_, config.includeDCDFile_, config.outDCDFile_);
    CSimulationBox simBox((CMolTwisterState*)state, stdOut, config);

    COut::setOutputFile(fopen(config.outInfoFile_.data(), "w"));
    mdLoop.runOptimization(simBox, config.gradientDescentAccuracy_, config.gradientDescentMaxSteps_, config.outputStride_);
    COut::closeOutputFile();
}

END_CUDA_COMPATIBLE()
