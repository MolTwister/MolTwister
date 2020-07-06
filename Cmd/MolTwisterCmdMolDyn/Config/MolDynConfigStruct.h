#pragma once
#include <string>

struct SMolDynConfigStruct
{
    enum Ensemble { ensembleNVT=0, ensembleNPT=1 };

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
};
