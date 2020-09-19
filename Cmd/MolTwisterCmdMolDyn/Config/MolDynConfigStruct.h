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
    std::string outPDistrFile_;
    std::string outVDistrFile_;
    bool includeXYZFile_;
    bool includeDCDFile_;
    bool includePDistrFile_;
    bool includeVDistrFile_;
    double maxPDistrOutput_;
    double maxVDistrOutput_;
    bool verboseOutput_;
    double scale12Interactions_;
    double scale13Interactions_;
    double scale14Interactions_;
    double scaleAbove14BondedInteractions_;
};
