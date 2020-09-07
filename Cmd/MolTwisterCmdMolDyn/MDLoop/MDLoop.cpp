#include "MDLoop.h"
#include "Printf.h"
#include "../Integrators/Constants.h"
#include "../../../Utilities/DCDFile.h"

BEGIN_CUDA_COMPATIBLE()

double CFctT::G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta)
{
    const double    g = double(SB_->dim_ * SB_->N_);
    const double    g_tilde = 1.0;
    double          ret = 0.0;
    
    
    if(j == 0)
    {
        double Sum_p2_m = 0.0;
        for(int k=0; k<SB_->N_; k++)
        {
            C3DVector p = SB_->particles_[k].p_;
            double p2 = p*p;
            Sum_p2_m+= (p2 / SB_->particles_[k].m_);
        }
        ret = Sum_p2_m - (g / beta);
    }
    else
    {
        ret = (p_eta[j-1]*p_eta[j-1] / Q[j-1]) - (g_tilde / beta);
    }
    
    return ret;
}

void CFctT::scaleMomentum(double coeff)
{
    for(int k=0; k<SB_->N_; k++)
    {
        SB_->particles_[k].p_*= coeff;
    }
}


double CFctP::G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta)
{
    const double g = 1.0;
    double ret = 0.0;
    
    if(j == 0)
    {
        ret = VV_->p_eps_*VV_->p_eps_ / VV_->W_ - (g / beta);
    }
    else
    {
        ret = (p_eta[j-1]*p_eta[j-1] / Q[j-1]) - (g / beta);
    }
    
    return ret;
}

void CFctP::scaleMomentum(double coeff)
{
    VV_->p_eps_*= coeff;
}



CMDLoop::CMDLoop(bool includeXYZFile, std::string fileNameXYZ, bool includeDCDFile, std::string fileNameDCD)
{
    includeXYZFile_ = includeXYZFile;
    includeDCDFile_ = includeDCDFile;
    fileNameXYZ_ = fileNameXYZ;
    fileNameDCD_ = fileNameDCD;

    const double maxV = 0.3; // Max v[AA/fs]in distr.
    maxP_ = maxV / Conv_P; // m*v_max, m = 1g/mol
}

void CMDLoop::runSimulation(CSimulationBox& simBox, int NStep, int outputEvery)
{
    CFctT fctT(&simBox);
    CFctP fctP(&simBox.velVerlet_);
    const int numPDistrBins = 150;
    const int equilibSteps = 1000;
    std::vector<int> momentumDistr[4];
    std::vector<int> volumeDistr;
    C3DVector boxSizeOut;
    mthost_vector<CMDFFMatrices::CForces>  F;

    srand((unsigned int)time(nullptr));
    
    resizeDistrArrays(momentumDistr, volumeDistr, numPDistrBins, 4);
    calcInitialForces(simBox, F);
    printHeading(simBox);
    
    for(int t=0; t<NStep; t++)
    {
        simBox.NHPPropagator(fctP);
        simBox.NHTPropagator(fctT);
        simBox.velVerPropagator(F, boxSizeOut);
        simBox.NHTPropagator(fctT);
        simBox.NHPPropagator(fctP);
        simBox.pbcWrap();

        updateOutput(t, equilibSteps, outputEvery, simBox, F, momentumDistr, volumeDistr, boxSizeOut);
    }
    
    finalizeOutput(simBox, momentumDistr, volumeDistr);
}

void CMDLoop::calcInitialForces(CSimulationBox& simBox, mthost_vector<CMDFFMatrices::CForces>& F)
{
    F = simBox.calcParticleForces();
}

void CMDLoop::printHeading(CSimulationBox& simBox)
{
    COut::printf("\r\n");
    COut::printf("\t----------------------------\r\n");
    COut::printf("\tMolecular dynamics simulation\r\n");
    COut::printf("\t----------------------------\r\n");
    COut::printf("\t Ensemble = %s\r\n", (simBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) ? "NPT" : ((simBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) ? "NVT" : "NVE"));
    if(simBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT || simBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) COut::printf("\t Temperature, T = %g K\r\n", simBox.NH_T_.T_*Conv_T);
    if(simBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) COut::printf("\t Pressure, P = %g atm\r\n", simBox.velVerlet_.P_*Conv_press);
    if(simBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT || simBox.getEnsemble() == SMolDynConfigStruct::ensembleNVT) COut::printf("\t NH T relax, tau = %g fs\r\n", simBox.NH_T_.tau_*Conv_t);
    if(simBox.getEnsemble() == SMolDynConfigStruct::ensembleNPT) COut::printf("\t NH P relax, tau = %g fs\r\n", simBox.NH_P_.tau_*Conv_t);
    COut::printf("\t Timestep, dt = %g fs\r\n", simBox.dt_*Conv_t);
    COut::printf("\t Particles, N = %i\r\n", simBox.N_);
    COut::printf("\t NH RESPA, n = %i\r\n", simBox.NH_T_.n_);
    COut::printf("\t NH Chain length, M = %i\r\n", simBox.NH_T_.M_);
    COut::printf("\t Init. box length, Lx = %g, Ly = %g, Lz = %g AA\r\n", simBox.getLmaxX(), simBox.getLmaxY(), simBox.getLmaxZ());
    COut::printf("\t Dimension, d = %i\r\n", simBox.dim_);
    COut::printf("\t----------------------------\r\n\r\n");
    COut::printf("\t%-15s%-15s%-15s%-20s\r\n", "Timestep", "Temp[K]", "Press[atm]", "Vol[AA^3]");
}

void CMDLoop::appendToXYZFile(mthost_vector<CParticle3D>& particles, int t, CSimulationBox& simBox)
{
    FILE* filePtr = fopen(fileNameXYZ_.data(), "a+");
    
    if(!filePtr) return;
    
    fprintf(filePtr, "%i\r\nFrame %u\r\n", (int)particles.size(), t);
    for(int k=0; k<(int)particles.size(); k++)
    {
        double x = particles[k].x_.x_;
        double y = particles[k].x_.y_;
        double z = particles[k].x_.z_;
        
        fprintf(filePtr, "%-15s% -15g%-15g%-15g\r\n", simBox.getAtomType(k).data(), x, y, z);
    }
    
    fclose(filePtr);
}

void CMDLoop::appendToDCDFile(mthost_vector<CParticle3D>& particles, CSimulationBox& simBox, const C3DVector& boxSize)
{
    FILE* file = fopen(fileNameDCD_.data(), "r");

    // Check if DCD file exists. If not, create it and add an appropriate header
    if(!file)
    {
        file = fopen(fileNameDCD_.data(), "w");
        if(file)
        {
            int numTimeSteps = simBox.getNumTimeSteps();
            int outputStride = simBox.getOutputStride();

            CDCDFile::CMainHeader mainHeader;
            mainHeader.ID_ = "CORD";
            mainHeader.nSets_ = numTimeSteps / outputStride;
            mainHeader.initStep_ = 0;
            mainHeader.wrtFreq_ = outputStride;
            mainHeader.timeStep_ = (float)simBox.getTimeStep();
            mainHeader.descriptA_ = "Written by MolTwister";
            mainHeader.descriptB_ = "---";
            mainHeader.nAtoms_ = (int)particles.size();

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
    record.init((int)particles.size(), true);

    for(size_t i=0; i<particles.size(); i++)
    {
        C3DVector r = particles[i].x_;
        record.setPos((int)i, r.x_, r.y_, r.z_);
    }
    record.write(file, numBytesWritten);

    fclose(file);
}

void CMDLoop::resizeDistrArrays(std::vector<int>* momentumDistr, std::vector<int>& volumeDistr, int size, int NArrays)
{
    for(int axis=0; axis<NArrays; axis++)
        momentumDistr[axis].resize(size, 0);
    
    volumeDistr.resize(size, 0);
}

void CMDLoop::appendToMomentumDistribution(CSimulationBox& simBox,
                                           std::vector<int>& momentumDistr,
                                           double maxP, int axis)
{
    double p;
    int iLen = (axis < 3) ? (int)simBox.particles_.size() : 1;
    
    for(int k=0; k<iLen; k++)
    {
        if(axis == 0)      p = simBox.particles_[k].p_.x_;
        else if(axis == 1) p = simBox.particles_[k].p_.y_;
        else if(axis == 2) p = simBox.particles_[k].p_.z_;
        else if(axis == 3) p = simBox.NH_T_.p_eta_[0];
        else                p = 0.0;
        
        int i = int(((p+maxP) / (2.0*maxP)) * double(momentumDistr.size()));
        if((i >=0) && (i < (int)momentumDistr.size())) momentumDistr[i]++;
    }
}

void CMDLoop::appendToVolumeDistribution(double V, std::vector<int>& volumeDistr, double maxV)
{
    int i = int((V / maxV) * double(volumeDistr.size()));
    if(i < (int)volumeDistr.size()) volumeDistr[i]++;
}

void CMDLoop::storeMomentumDistribution(std::string fileName,
                                        std::vector<int>& momentumDistr,
                                        double maxP, int axis)
{
    double convFact = (axis < 3) ? Conv_P : Conv_Peta;
    FILE* file = fopen(fileName.data(), "w");
    int N = (int)momentumDistr.size();
    
    if(!file) return;

    std::string szUnit = (axis < 3) ? "p[gAA/(mol*fs)]" : "p[kJ*fs/mol]";
    fprintf(file, "%-20s%-20s\r\n", szUnit.data(), "P");
    for(int i=0; i<N; i++)
    {
        double dN = double(N);
        double di = double(i);
        double binCent = ((maxP / dN)*(2.0*di + 1.0) - maxP) * convFact;
        fprintf(file, "%-20g%-20i\r\n", binCent, momentumDistr[i]);
    }
    
    fclose(file);
}

void CMDLoop::storeVolumeDistribution(std::string fileName, std::vector<int>& volumeDistr, double maxV)
{
    FILE* filePtr = fopen(fileName.data(), "w");
    int N = (int)volumeDistr.size();
    
    if(!filePtr) return;

    fprintf(filePtr, "%-20s%-20s\r\n", "V[AA^3]", "P");
    for(int i=0; i<N; i++)
    {
        double dBinCent = (maxV / (2.0*double(N))) * double(2*i + 1);
        fprintf(filePtr, "%-20g%-20i\r\n", dBinCent, volumeDistr[i]);
    }

    fclose(filePtr);
}

void CMDLoop::updateOutput(int t, int equilibSteps, int outputEvery, CSimulationBox& simBox,
                           const mthost_vector<CMDFFMatrices::CForces>& F, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr,
                           const C3DVector& boxSize)
{
    // Perform calculations on MD trajectories and output data
    if(t > equilibSteps)
    {
        for(int axis=0; axis<3; axis++)
        {
            appendToMomentumDistribution(simBox, momentumDistr[axis], maxP_, axis);
        }
        
        double maxPeta  = 2.0*simBox.NH_T_.Q_[0]/simBox.NH_T_.tau_;
        appendToMomentumDistribution(simBox, momentumDistr[3], maxPeta, 3);

        double maxV = 1000.0E3;
        appendToVolumeDistribution(simBox.calcV(), volumeDistr, maxV);
    }
    
    if((t % outputEvery) == 0)
    {
        if(includeXYZFile_) appendToXYZFile(simBox.particles_, t, simBox);
        if(includeDCDFile_) appendToDCDFile(simBox.particles_, simBox, boxSize);

        double T = simBox.calcTemp() * Conv_T;
        COut::printf("\t%-15i%-15g%-15g%-20g\r\n", t, T,
               simBox.calcPress(F) * Conv_press, simBox.calcV());
    }
}

void CMDLoop::finalizeOutput(CSimulationBox& simBox, std::vector<int>* momentumDistr, std::vector<int>& volumeDistr)
{
    // Output final momentum distribution
    double maxPeta = 2.0*simBox.NH_T_.Q_[0]/simBox.NH_T_.tau_;
    double maxV = 1000.0E3;
    storeMomentumDistribution("gas_distr_x.dat", momentumDistr[0], maxP_, 0);
    storeMomentumDistribution("gas_distr_y.dat", momentumDistr[1], maxP_, 1);
    storeMomentumDistribution("gas_distr_z.dat", momentumDistr[2], maxP_, 2);
    storeMomentumDistribution("gas_distr_p_eta0.dat", momentumDistr[3], maxPeta, 3);
    storeVolumeDistribution("gas_distr_v.dat", volumeDistr, maxV);
}

END_CUDA_COMPATIBLE()
