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
    CFctT(CSimulationBox* pSB) : CFct() { m_pSB = pSB; }
    
public:
    double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta);
    void ScaleMomentum(double coeff);
    
private:
    CSimulationBox* m_pSB;
};

class CFctP : public CFct
{
public:
    CFctP(CVelVerlet* pVV) : CFct() { m_pVV = pVV; }
    
public:
    double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta);
    void ScaleMomentum(double coeff);
    
private:
    CVelVerlet* m_pVV;
};

class CMDLoop
{
public:
    CMDLoop();

public:
    void RunSimulation(CSimulationBox& SimBox, int iNStep, int iOutputEvery);
    
private:
    void CalcInitialForces(CSimulationBox& SimBox, mthost_vector<CMDFFMatrices::CForces>& F);
    void NegMomHalfWay(int t, int iNStep, CSimulationBox& SimBox);
    void PrintHeading(CSimulationBox& SimBox);
    void AppendToXYZFile(mthost_vector<CParticle3D>& aParticles, int t);
    void ResizeDistrArrays(std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr, int iSize, int iNArrays);
    void AppendToMomentumDistribution(CSimulationBox& SimBox, std::vector<int>& aMomentumDistr, double dMaxP, int iAxis);
    void AppendToVolumeDistribution(double V, std::vector<int>& aVolumeDistr, double dMaxV);
    void StoreMomentumDistribution(std::string szFileName, std::vector<int>& aMomentumDistr, double dMaxP, int iAxis);
    void StoreVolumeDistribution(std::string szFileName, std::vector<int>& aVolumeDistr, double dMaxV);
    void UpdateOutput(int t, int iEquilibSteps, int iOutputEvery, CSimulationBox& SimBox, const mthost_vector<CMDFFMatrices::CForces>& F, std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr);
    void FinalizeOutput(CSimulationBox& SimBox, std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr);
    
private:
    double m_dMaxP;
};

END_CUDA_COMPATIBLE()
