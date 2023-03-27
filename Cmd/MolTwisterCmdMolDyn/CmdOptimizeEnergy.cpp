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

#include "CmdOptimizeEnergy.h"
#include "../../Utilities/ASCIIUtility.h"
#include "Simulator/MDSimulator.h"

#if INCLUDE_CUDA_COMMANDS == 1
namespace mtdev
{
    #include "Simulator/MDSimulator.h"
}
#endif

std::string CCmdOptimizeEnergy::getCmd()
{
    return "optimizeenergy";
}

std::vector<std::string> CCmdOptimizeEnergy::getCmdLineKeywords()
{
    return { "optimizeenergy", "cpu" };
}

std::vector<std::string> CCmdOptimizeEnergy::getCmdHelpLines()
{
    return {
        "optimizeenergy [cpu]"
    };
}

std::string CCmdOptimizeEnergy::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tThis command will run an energy optimization. If the 'cpu' keyword is included, the\r\n";
    text+= "\toptimization will be forced to be executed on the CPU. If not, an attempt will be made to run the\r\n";
    text+= "\toptimization on GPU. If this does not succeed (e.g., if the software is not compiled to run on GPU)\r\n";
    text+= "\tthe optimization will be executed on CPU.";

    return text;
}

std::string CCmdOptimizeEnergy::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;

    // Make sure our MD system is centralized
    int atomCount = (int)state_->atoms_.size();
    C3DRect boundingBox;
    std::vector<int> atomIndices(atomCount);
    for(int i=0; i<atomCount; i++) atomIndices[i] = i;
    boundingBox = state_->calcBoundingBox(state_->currentFrame_, atomIndices);
    C3DVector v0;
    for(int i=0; i<atomCount; i++)
    {
        CAtom* atomPtr = state_->atoms_[i].get();

        v0.x_ = (boundingBox.rLow_.x_ + boundingBox.rHigh_.x_) / 2.0;
        v0.y_ = (boundingBox.rLow_.y_ + boundingBox.rHigh_.y_) / 2.0;
        v0.z_ = (boundingBox.rLow_.z_ + boundingBox.rHigh_.z_) / 2.0;
        atomPtr->r_[state_->currentFrame_]+= v0;
    }

    // Determine if simulation should be executed on CPU or on GPU
    std::string text = CASCIIUtility::getArg(arguments, arg++);
    bool runOnGPU = false;
    bool compiledForGPU = false;
    #if INCLUDE_CUDA_COMMANDS == 1
    {
        compiledForGPU = true;
    }
    #endif

    if(text == "cpu")
    {
        runOnGPU = false;
        if(compiledForGPU)
        {
            fprintf(stdOut_, "\tMolTwister is compiled for GPU, but 'cpu' was explicitly specified. Hence, energy optimization will run on CPU!\r\n");
        }
        else
        {
            fprintf(stdOut_, "\tMolTwister is NOT compiled for GPU, but 'cpu' was explicitly specified. Hence, engergy optimization will run on CPU!\r\n");
        }
    }
    else
    {
        if(compiledForGPU)
        {
            fprintf(stdOut_, "\tMolTwister is compiled for GPU and will run energy optimization on GPU!\r\n");
            runOnGPU = true;
        }
        else
        {
            fprintf(stdOut_, "\tMolTwister is NOT compiled for GPU and will therefore run energy optimization on CPU!\r\n");
            runOnGPU = false;
        }
    }

    // Run simulation
    if(runOnGPU)
    {
        #if INCLUDE_CUDA_COMMANDS == 1
        {
            CSerializer stateContent;
            state_->serialize(stateContent, true);
            mtdev::CMDSimulator::optimize(molDynConfig_->cfg_, stdOut_, stateContent);
        }
        #else
        {
            lastError_ = "\tMolTwister is NOT compiled for GPU, but the simulation attemptet to run on GPU!\r\n";
            return lastError_;
        }
        #endif
    }
    else
    {
        CMDSimulator::optimize(molDynConfig_->cfg_, stdOut_, (CMolTwisterState*)state_);
    }

    return lastError_;
}

