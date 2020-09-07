#pragma once
#include <string>
#include <vector>
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
public:
    CMDLoop() = delete;
    CMDLoop(bool includeXYZFile, std::string fileNameXYZ, bool includeDCDFile, std::string fileNameDCD);

public:
    void runSimulation(CSimulationBox& simBox, int NStep, int outputEvery);
    
private:
    void calcInitialForces(CSimulationBox& simBox, mthost_vector<CMDFFMatrices::CForces>& F);
    void printHeading(CSimulationBox& simBox);
    void appendToXYZFile(mthost_vector<CParticle3D>& particles, int t, CSimulationBox& simBox);
    void appendToDCDFile(mthost_vector<CParticle3D>& particles, CSimulationBox& simBox, const C3DVector& boxSize);
    void resizeDistrArrays(std::vector<int>* momentumDistr, std::vector<int>& volumeDistr, int size, int NArrays);
    void appendToMomentumDistribution(CSimulationBox& simBox, std::vector<int>& momentumDistr, double maxP, int axis);
    void appendToVolumeDistribution(double V, std::vector<int>& volumeDistr, double maxV);
    void storeMomentumDistribution(std::string fileName, std::vector<int>& momentumDistr, double maxP, int axis);
    void storeVolumeDistribution(std::string fileName, std::vector<int>& volumeDistr, double maxV);
    void updateOutput(int t, int equilibSteps, int outputEvery, CSimulationBox& simBox, const mthost_vector<CMDFFMatrices::CForces>& F, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr, const C3DVector& boxSize);
    void finalizeOutput(CSimulationBox& simBox, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr);
    
private:
    double maxP_;
    bool includeXYZFile_;
    bool includeDCDFile_;
    std::string fileNameXYZ_;
    std::string fileNameDCD_;
};

END_CUDA_COMPATIBLE()
