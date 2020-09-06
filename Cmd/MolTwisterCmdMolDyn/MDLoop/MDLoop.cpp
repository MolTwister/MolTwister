#include "MDLoop.h"
#include "Printf.h"
#include "../Integrators/Constants.h"
#include "../../../Utilities/DCDFile.h"

BEGIN_CUDA_COMPATIBLE()

double CFctT::G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta)
{
    const double    g = double(m_pSB->dim * m_pSB->N);
    const double    g_tilde = 1.0;
    double          ret = 0.0;
    
    
    if(j == 0)
    {
        double Sum_p2_m = 0.0;
        for(int k=0; k<m_pSB->N; k++)
        {
            C3DVector p = m_pSB->aParticles[k].p;
            double p2 = p*p;
            Sum_p2_m+= (p2 / m_pSB->aParticles[k].m);
        }
        ret = Sum_p2_m - (g / beta);
    }
    else
    {
        ret = (p_eta[j-1]*p_eta[j-1] / Q[j-1]) - (g_tilde / beta);
    }
    
    return ret;
}

void CFctT::ScaleMomentum(double coeff)
{
    for(int k=0; k<m_pSB->N; k++)
    {
        m_pSB->aParticles[k].p*= coeff;
    }
}


double CFctP::G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta)
{
    const double    g = 1.0;
    double          ret = 0.0;
    
    
    if(j == 0)
    {
        ret = m_pVV->p_eps*m_pVV->p_eps / m_pVV->W - (g / beta);
    }
    else
    {
        ret = (p_eta[j-1]*p_eta[j-1] / Q[j-1]) - (g / beta);
    }
    
    return ret;
}

void CFctP::ScaleMomentum(double coeff)
{
    m_pVV->p_eps*= coeff;
}



CMDLoop::CMDLoop(std::string fileNameXYZ)
{
    fileNameXYZ_ = fileNameXYZ;
    fileNameDCD_ = "traj.dcd";

    const double dMaxV = 0.3; // Max v[Å/fs]in distr.
    m_dMaxP = dMaxV / Conv_P; // m*v_max, m = 1g/mol
}

void CMDLoop::RunSimulation(CSimulationBox& SimBox, int iNStep, int iOutputEvery)
{
    CFctT                     FctT(&SimBox);
    CFctP                     FctP(&SimBox.VelVerlet);
    const int                 iNumPDistrBins = 150;
    const int                 iEquilibSteps = 1000;
    std::vector<int>          aMomentumDistr[4];
    std::vector<int>          aVolumeDistr;
    C3DVector                 boxSizeOut;
    mthost_vector<CMDFFMatrices::CForces>  F;
    
    
    srand((unsigned int)time(nullptr));
    
    ResizeDistrArrays(aMomentumDistr, aVolumeDistr, iNumPDistrBins, 4);
    CalcInitialForces(SimBox, F);
    PrintHeading(SimBox);
    
    for(int t=0; t<iNStep; t++)
    {
        SimBox.NHPPropagator(FctP);
        SimBox.NHTPropagator(FctT);
        SimBox.VelVerPropagator(F, boxSizeOut);
        SimBox.NHTPropagator(FctT);
        SimBox.NHPPropagator(FctP);
        SimBox.PBCWrap();

        UpdateOutput(t, iEquilibSteps, iOutputEvery, SimBox, F, aMomentumDistr, aVolumeDistr, boxSizeOut);
    }
    
    FinalizeOutput(SimBox, aMomentumDistr, aVolumeDistr);
}

void CMDLoop::CalcInitialForces(CSimulationBox& SimBox, mthost_vector<CMDFFMatrices::CForces>& F)
{
    F = SimBox.CalcParticleForces();
}

void CMDLoop::PrintHeading(CSimulationBox& SimBox)
{
    COut::Printf("\r\n");
    COut::Printf("\t----------------------------\r\n");
    COut::Printf("\tMolecular dynamics simulation\r\n");
    COut::Printf("\t----------------------------\r\n");
    COut::Printf("\t Ensemble = %s\r\n", (SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) ? "NPT" : ((SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) ? "NVT" : "NVE"));
    if(SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT || SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) COut::Printf("\t Temperature, T = %g K\r\n", SimBox.NH_T.T*Conv_T);
    if(SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) COut::Printf("\t Pressure, P = %g atm\r\n", SimBox.VelVerlet.P*Conv_press);
    if(SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT || SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) COut::Printf("\t NH T relax, tau = %g fs\r\n", SimBox.NH_T.tau*Conv_t);
    if(SimBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) COut::Printf("\t NH P relax, tau = %g fs\r\n", SimBox.NH_P.tau*Conv_t);
    COut::Printf("\t Timestep, dt = %g fs\r\n", SimBox.dt*Conv_t);
    COut::Printf("\t Particles, N = %i\r\n", SimBox.N);
    COut::Printf("\t NH RESPA, n = %i\r\n", SimBox.NH_T.n);
    COut::Printf("\t NH Chain length, M = %i\r\n", SimBox.NH_T.M);
    COut::Printf("\t Init. box length, Lx = %g, Ly = %g, Lz = %g Å\r\n", SimBox.getLmaxX(), SimBox.getLmaxY(), SimBox.getLmaxZ());
    COut::Printf("\t Dimension, d = %i\r\n", SimBox.dim);
    COut::Printf("\t----------------------------\r\n\r\n");
    COut::Printf("\t%-15s%-15s%-15s%-20s\r\n", "Timestep", "Temp[K]", "Press[atm]", "Vol[Å^3]");
}

void CMDLoop::AppendToXYZFile(mthost_vector<CParticle3D>& aParticles, int t, CSimulationBox& SimBox)
{
    FILE* pFile = fopen(fileNameXYZ_.data(), "a+");
    
    if(!pFile) return;
    
    fprintf(pFile, "%i\r\nFrame %u\r\n", (int)aParticles.size(), t);
    for(int k=0; k<(int)aParticles.size(); k++)
    {
        double x = aParticles[k].x.x_;
        double y = aParticles[k].x.y_;
        double z = aParticles[k].x.z_;
        
        fprintf(pFile, "%-15s% -15g%-15g%-15g\r\n", SimBox.getAtomType(k).data(), x, y, z);
    }
    
    fclose(pFile);
}

void CMDLoop::appendToDCDFile(mthost_vector<CParticle3D>& aParticles, CSimulationBox& SimBox, const C3DVector& boxSize)
{
    FILE* file = fopen(fileNameDCD_.data(), "r");

    // Check if DCD file exists. If not, create it and add an appropriate header
    if(!file)
    {
        file = fopen(fileNameDCD_.data(), "w");
        if(file)
        {
            int numTimeSteps = SimBox.getNumTimeSteps();
            int outputStride = SimBox.getOutputStride();

            CDCDFile::CMainHeader mainHeader;
            mainHeader.ID_ = "CORD";
            mainHeader.nSets_ = numTimeSteps / outputStride;
            mainHeader.initStep_ = 0;
            mainHeader.wrtFreq_ = outputStride;
            mainHeader.timeStep_ = (float)SimBox.getTimeStep();
            mainHeader.descriptA_ = "Written by MolTwister";
            mainHeader.descriptB_ = "---";
            mainHeader.nAtoms_ = (int)aParticles.size();

            int numBytesWritten = 0;
            mainHeader.write(file, numBytesWritten);
        }

        fclose(file);
    }
    else
    {
        fclose(file);
    }

    // Open the DCD file and place the file pointer at the end of the file
    file = fopen(fileNameDCD_.data(), "a+");
    if(!file) return;

    // Write a DCD record, starting at the end of the file
    CDCDFile::CRecordHeader recordHeader;
    recordHeader.boxX_ = boxSize.x_;
    recordHeader.boxY_ = boxSize.y_;
    recordHeader.boxZ_ = boxSize.z_;

    CDCDFile::CRecord record;
    record.setRecordHeader(recordHeader);

    int numBytesWritten = 0;
    record.init((int)aParticles.size(), true);

    for(size_t i=0; i<aParticles.size(); i++)
    {
        C3DVector r = aParticles[i].x;
        record.setPos((int)i, r.x_, r.y_, r.z_);
    }
    record.write(file, numBytesWritten);

    fclose(file);
}

void CMDLoop::ResizeDistrArrays(std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr, int iSize, int iNArrays)
{
    for(int iAxis=0; iAxis<iNArrays; iAxis++)
        aMomentumDistr[iAxis].resize(iSize, 0);
    
    aVolumeDistr.resize(iSize, 0);
}

void CMDLoop::AppendToMomentumDistribution(CSimulationBox& SimBox,
                                           std::vector<int>& aMomentumDistr,
                                           double dMaxP, int iAxis)
{
    double  p;
    int     iLen = (iAxis < 3) ? (int)SimBox.aParticles.size() : 1;
    
    for(int k=0; k<iLen; k++)
    {
        if(iAxis == 0)      p = SimBox.aParticles[k].p.x_;
        else if(iAxis == 1) p = SimBox.aParticles[k].p.y_;
        else if(iAxis == 2) p = SimBox.aParticles[k].p.z_;
        else if(iAxis == 3) p = SimBox.NH_T.p_eta[0];
        else                p = 0.0;
        
        int i = int(((p+dMaxP) / (2.0*dMaxP)) * double(aMomentumDistr.size()));
        if((i >=0) && (i < (int)aMomentumDistr.size())) aMomentumDistr[i]++;
    }
}

void CMDLoop::AppendToVolumeDistribution(double V, std::vector<int>& aVolumeDistr, double dMaxV)
{
    int i = int((V / dMaxV) * double(aVolumeDistr.size()));
    if(i < (int)aVolumeDistr.size()) aVolumeDistr[i]++;
}

void CMDLoop::StoreMomentumDistribution(std::string szFileName,
                                        std::vector<int>& aMomentumDistr,
                                        double dMaxP, int iAxis)
{
    double      dConvFact = (iAxis < 3) ? Conv_P : Conv_Peta;
    FILE*       pFile = fopen(szFileName.data(), "w");
    int         iN = (int)aMomentumDistr.size();
    
    if(!pFile) return;

    std::string szUnit = (iAxis < 3) ? "p[gÅ/(mol*fs)]" : "p[kJ*fs/mol]";
    fprintf(pFile, "%-20s%-20s\r\n", szUnit.data(), "P");
    for(int i=0; i<iN; i++)
    {
        double dN = double(iN);
        double di = double(i);
        double dBinCent = ((dMaxP / dN)*(2.0*di + 1.0) - dMaxP) * dConvFact;
        fprintf(pFile, "%-20g%-20i\r\n", dBinCent, aMomentumDistr[i]);
    }
    
    fclose(pFile);
}

void CMDLoop::StoreVolumeDistribution(std::string szFileName, std::vector<int>& aVolumeDistr, double dMaxV)
{
    FILE*       pFile = fopen(szFileName.data(), "w");
    int         iN = (int)aVolumeDistr.size();
    
    if(!pFile) return;

    fprintf(pFile, "%-20s%-20s\r\n", "V[AA^3]", "P");
    for(int i=0; i<iN; i++)
    {
        double dBinCent = (dMaxV / (2.0*double(iN))) * double(2*i + 1);
        fprintf(pFile, "%-20g%-20i\r\n", dBinCent, aVolumeDistr[i]);
    }

    fclose(pFile);
}

void CMDLoop::UpdateOutput(int t, int iEquilibSteps, int iOutputEvery, CSimulationBox& SimBox,
                           const mthost_vector<CMDFFMatrices::CForces>& F, std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr,
                           const C3DVector& boxSize)
{
    // Perform calculations on MD trajectories and output data
    if(t > iEquilibSteps)
    {
        for(int iAxis=0; iAxis<3; iAxis++)
        {
            AppendToMomentumDistribution(SimBox, aMomentumDistr[iAxis], m_dMaxP, iAxis);
        }
        
        double dMaxPeta  = 2.0*SimBox.NH_T.Q[0]/SimBox.NH_T.tau;
        AppendToMomentumDistribution(SimBox, aMomentumDistr[3], dMaxPeta, 3);

        double dMaxV = 1000.0E3;
        AppendToVolumeDistribution(SimBox.CalcV(), aVolumeDistr, dMaxV);
    }
    
    if((t % iOutputEvery) == 0)
    {
        // :TODO: Should have an option to choose either XYZ, DCD or both
        //        also remember to make the DCD name configurable
        AppendToXYZFile(SimBox.aParticles, t, SimBox);
        appendToDCDFile(SimBox.aParticles, SimBox, boxSize);
        double T = SimBox.CalcTemp() * Conv_T;
        COut::Printf("\t%-15i%-15g%-15g%-20g\r\n", t, T,
               SimBox.CalcPress(F) * Conv_press, SimBox.CalcV());
    }
}

void CMDLoop::FinalizeOutput(CSimulationBox& SimBox, std::vector<int>* aMomentumDistr, std::vector<int>& aVolumeDistr)
{
    // Output final momentum distribution
    double dMaxPeta = 2.0*SimBox.NH_T.Q[0]/SimBox.NH_T.tau;
    double dMaxV = 1000.0E3;
    StoreMomentumDistribution("gas_distr_x.dat", aMomentumDistr[0], m_dMaxP, 0);
    StoreMomentumDistribution("gas_distr_y.dat", aMomentumDistr[1], m_dMaxP, 1);
    StoreMomentumDistribution("gas_distr_z.dat", aMomentumDistr[2], m_dMaxP, 2);
    StoreMomentumDistribution("gas_distr_p_eta0.dat", aMomentumDistr[3], dMaxPeta, 3);
    StoreVolumeDistribution("gas_distr_v.dat", aVolumeDistr, dMaxV);
}

END_CUDA_COMPATIBLE()
