#ifndef __ThesisMDTests__MDLoop__
#define __ThesisMDTests__MDLoop__

#include "../SimulationBox/SimulationBox.h"
#include "../Integrators/VelVerlet.h"
#include <string>
#include <vector>

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
    void CalcInitialForces(CSimulationBox& SimBox, vector<CMDForces>& F);
    void NegMomHalfWay(int t, int iNStep, CSimulationBox& SimBox);
    void PrintHeading(CSimulationBox& SimBox);
    void AppendToXYZFile(vector<CParticle3D>& aParticles, int t);
    void ResizeDistrArrays(vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr, int iSize, int iNArrays);
    void AppendToMomentumDistribution(CSimulationBox& SimBox, std::vector<int>& aMomentumDistr, double dMaxP, int iAxis);
    void AppendToVolumeDistribution(double V, vector<int>& aVolumeDistr, double dMaxV);
    void StoreMomentumDistribution(std::string szFileName, std::vector<int>& aMomentumDistr, double dMaxP, int iAxis);
    void StoreVolumeDistribution(std::string szFileName, std::vector<int>& aVolumeDistr, double dMaxV);
    void UpdateOutput(int t, int iEquilibSteps, int iOutputEvery, CSimulationBox& SimBox, const vector<CMDForces>& F, vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr);
    void FinalizeOutput(CSimulationBox& SimBox, std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr);
    
private:
    double m_dMaxP;
};

#endif
