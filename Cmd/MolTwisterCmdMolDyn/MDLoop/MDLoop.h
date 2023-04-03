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
#include <vector>
#include <ctime>
#include <chrono>
#include "../SimulationBox/SimulationBox.h"
#include "../Integrators/VelVerlet.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFctT : public CFct
{
public:
    CFctT(CSimulationBox* SB) : CFct() { SB_ = SB; }
    virtual ~CFctT() {}
    
public:
    double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta);
    void scaleMomentum(double coeff);
    
private:
    CSimulationBox* SB_;
};

class CFctP : public CFct
{
public:
    CFctP(CVelVerlet* VV) : CFct() { VV_ = VV; }
    virtual ~CFctP() {}

public:
    double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta);
    void scaleMomentum(double coeff);
    
private:
    CVelVerlet* VV_;
};

class CMDLoop
{
private:
    enum EHeading { headingMD=0, headingOptimization };
public:
    CMDLoop() = delete;
    CMDLoop(bool includeXYZFile, std::string fileNameXYZ, bool includeDCDFile, std::string fileNameDCD, bool includeXTCFile, std::string fileNameXTC, const float& xtcPrecision);

public:
    void runSimulation(CSimulationBox& simBox, int NStep, int outputEvery);
    void runOptimization(CSimulationBox& simBox, double accuracy, int maxSteps, int outputEvery);
    void setPDistrOutput(bool includePDistrOutput, double maxPDistrOutput, std::string fileNamePDistr);
    void setVDistrOutput(bool includeVDistrOutput, double setVDistrOutput, std::string fileNameVDistr);
    
private:
    void calcInitialForces(CSimulationBox& simBox, mthost_vector<CMDFFMatrices::CForces>& F);
    void printHeading(CSimulationBox& simBox, EHeading heading);
    void appendToXYZFile(mthost_vector<CParticle3D>& particles, int t, CSimulationBox& simBox) const;
    void appendToDCDFile(mthost_vector<CParticle3D>& particles, CSimulationBox& simBox, const C3DVector& boxSize, const int& numTimeSteps, const int& outputStride) const;
    void appendToXTCFile(mthost_vector<CParticle3D>& particles, CSimulationBox& simBox, const C3DVector& boxSize, const int& t, const int& numTimeSteps, const int& outputStride, const float& precision) const;
    void resizeDistrArrays(std::vector<int>* momentumDistr, std::vector<int>& volumeDistr, int size, int NArrays);
    void appendToMomentumDistribution(CSimulationBox& simBox, std::vector<int>& momentumDistr, double maxP, int axis);
    void appendToVolumeDistribution(double V, std::vector<int>& volumeDistr, double maxV);
    void storeMomentumDistribution(std::string fileName, std::vector<int>& momentumDistr, double maxP, int axis);
    void storeVolumeDistribution(std::string fileName, std::vector<int>& volumeDistr, double maxV);
    void updateOutput(int t, int totalSteps, int equilibSteps, int outputEvery, CSimulationBox& simBox, const mthost_vector<CMDFFMatrices::CForces>& F, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr, const C3DVector& boxSize);
    void finalizeOutput(CSimulationBox& simBox, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr);
    double calcMaxPetaDistrOutput(CSimulationBox& simBox);
    void startTimer();
    void endTimer();
    double readTimerAverage();

private:
    double maxPDistrOutput_;
    double maxVDistrOutput_;
    bool includeXYZFile_;
    bool includeDCDFile_;
    bool includeXTCFile_;
    bool includePDistrOutput_;
    bool includeVDistrOutput_;
    float xtcPrecision_;
    std::string fileNameXYZ_;
    std::string fileNameDCD_;
    std::string fileNameXTC_;
    std::string fileNamePDistr_;
    std::string fileNameVDistr_;
    std::chrono::steady_clock::time_point lastStartingTimeStep_;
    std::vector<std::chrono::duration<double>> timeStepTimePeriods_;
};

END_CUDA_COMPATIBLE()
