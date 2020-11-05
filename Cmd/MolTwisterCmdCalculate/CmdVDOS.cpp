#include "CmdVDOS.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../../Utilities/FFT1D.h"
#include "../Tools/ProgressBar.h"

std::string CCmdVDOS::getCmd()
{
    return "vdos";
}

std::vector<std::string> CCmdVDOS::getCmdLineKeywords()
{
    return { "vdos", "name", "sel" };
}

std::vector<std::string> CCmdVDOS::getCmdHelpLines()
{
    return {
                "vdos <DCD filename> <frame from> <num. bits> <time step (fs)> name <atom IDs (comma sep., no space)>",
                "vdos <DCD filename> <frame from> <num. bits> <time step (fs)> sel"
           };
}

std::string CCmdVDOS::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the vibrational density of states (VDOS) for the DCD file, <DCD filename>,\r\n";
    text+= "\tbetween the start frame, <frame from>, to the frame <frame from>+2^<num. bits>. The time\r\n";
    text+= "\tstep, <time step>, given in units of femtoseconds (fs), is used to obtain numerically\r\n";
    text+= "\tcalculated velocities that are needed within the calculations. VDOS is only calculated for a\r\n";
    text+= "\tgiven selection of atoms. Either, based on the atom 'name', where a list, <atom IDs>,\r\n";
    text+= "\tis supplied (e.g., H, O, C7), or based on the visual selection of atoms, achieved through\r\n";
    text+= "\tthe 'sel' keyword.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. omegaomega[rad/fs] vdos\r\n";
    text+= "\t2. <angular frequency> <VDOS value (probability)>\r\n";
    text+= "\t3. <angular frequency> <VDOS value (probability)>\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\tN+1. <angular frequency> <VDOS value (probability)>\r\n";
    text+= "\twhere N is the length of the VDOS function.";

    return text;
}

std::string CCmdVDOS::execute(std::vector<std::string> arguments)
{
    /////////////////////////////////////////////////////////////////////////////////
    // Calculate the vibrational density of states (VDOS) by the Fourier components
    // of the VACF (which can be found through FFT of the velocity trajectories).
    // Ref. [Phys. Rev. B, 46, 12019 (1992)], Ref. [Phys. Rev., 188, 1407 (1969)]
    // and Ref. [Statistical Mechanics: Theory and Molecular Simulation, Tuckerman]
    /////////////////////////////////////////////////////////////////////////////////
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToInclude;
    CDCDFile dcdFile;
    std::string text;
    double timeStep;
    int frameFrom, frameTo;


    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    frameFrom = atoi(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    frameTo = frameFrom + (int)pow(2.0, atoi(text.data()));
    text = CASCIIUtility::getArg(arguments, arg++);
    timeStep = atof(text.data());
    if(timeStep < 0.0)
    {
        lastError_ = "Error: timestep cannot be zero!";
        return lastError_;
    }


    // Find indices to loop over
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "name")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        CASCIIUtility::removeWhiteSpace(text);
        atomsToInclude = CASCIIUtility::getWords(text, ",");

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            text = state_->atoms_[i]->getID();
            for(int j=0; j<atomsToInclude.size(); j++)
            {
                if(text == atomsToInclude[j])
                {
                    atomIndicesToInclude.emplace_back(i);
                }
            }
        }
    }
    else if(text == "sel")
    {
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->isSelected())
            {
                atomIndicesToInclude.emplace_back(i);
            }
        }
    }
    else
    {
        lastError_ = "Error: expected 'name' or 'sel'!";
        return lastError_;
    }


    // Pre-read all selected atom trajectories (that need to be fftLen+1, while velocity traj. becomes fftLen long)
    int fftLen = frameTo - frameFrom;
    int numAtoms = (int)atomIndicesToInclude.size();
    int numFrames = dcdFile.getNumRecords();

    std::vector<std::vector<std::vector<CFFT1D::CCplx>>> atomTrajectory;
    atomTrajectory.resize(3);
    for(int n=0; n<3; n++)
    {
        atomTrajectory[n].resize(numAtoms);
        for(int i=0; i<numAtoms; i++)
        {
            atomTrajectory[n][i].resize(fftLen+1);
        }
    }

    pb. beginProgress("Reading atomic trajectories");
    for(int t=frameFrom; t<(frameTo+1); t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        C3DVector v;
        dcdFile.gotoRecord(t);
        for(int i=0; i<numAtoms; i++)
        {
            v = dcdFile.getCoordinate(atomIndicesToInclude[i]);
            atomTrajectory[0][i][t-frameFrom].re_ = v.x_;
            atomTrajectory[1][i][t-frameFrom].re_ = v.y_;
            atomTrajectory[2][i][t-frameFrom].re_ = v.z_;
        }

        pb.updateProgress(t-frameFrom, fftLen+1);
    }
    pb.endProgress();


    // Convert positional trajectories to velocity trajectories
    pb. beginProgress("Reading atomic trajectories");
    size_t trajAtomCount = atomTrajectory[0].size();
    if(trajAtomCount != atomTrajectory[1].size())
    {
        lastError_ = "The atom count of the y-trajectories does not match that of x";
        return lastError_;
    }
    if(trajAtomCount != atomTrajectory[2].size())
    {
        lastError_ = "The atom count of the z-trajectories does not match that of x";
        return lastError_;
    }

    if(trajAtomCount == 0)
    {
        lastError_ = "The atom count was zero";
        return lastError_;
    }

    size_t trajSize = atomTrajectory[0][0].size();
    if(trajSize != atomTrajectory[1][0].size())
    {
        lastError_ = "The trajectory size of the y-trajectories does not match that of x";
        return lastError_;
    }
    if(trajSize != atomTrajectory[2][0].size())
    {
        lastError_ = "The trajectory size of the z-trajectories does not match that of x";
        return lastError_;
    }

    std::vector<std::vector<CFFT1D::CCplx>> atomVelTrajectory(trajAtomCount);
    for(size_t i=0; i<atomVelTrajectory.size(); i++)
    {
        atomVelTrajectory[i].resize(fftLen);
    }

    for(size_t t=1; t<trajSize; t++)
    {
        for(size_t i=0; i<trajAtomCount; i++)
        {
            double x1 = atomTrajectory[0][i][t-1].re_;
            double x2 = atomTrajectory[0][i][t].re_;
            double vx = (x2 - x1) / timeStep; // AA / fs

            double y1 = atomTrajectory[1][i][t-1].re_;
            double y2 = atomTrajectory[1][i][t].re_;
            double vy = (y2 - y1) / timeStep; // AA / fs

            double z1 = atomTrajectory[2][i][t-1].re_;
            double z2 = atomTrajectory[2][i][t].re_;
            double vz = (z2 - z1) / timeStep; // AA / fs

            double kern = vx*vx + vy*vy + vz*vz;
            kern = (kern >= 0.0) ? kern : 0.0;
            double v = sqrt(kern);

            atomVelTrajectory[i][t-1].re_ = v;
        }

        pb.updateProgress((int)t, (int)trajSize);
    }

    pb.endProgress();


    // Loop through each particle, then FFT each of the particle velocity trajectories
    std::vector<double> avg_v_tilde2;
    avg_v_tilde2.resize(fftLen, 0.0);
    pb.beginProgress("Calculating vibrational density of states");
    for(int i=0; i<(int)trajAtomCount; i++)
    {
        // Fourier transform velocity trajectory of atom i
        CFFT1D fft;
        std::shared_ptr<std::vector<CFFT1D::CCplx>> v_tilde;
        v_tilde = fft.fft1D(atomVelTrajectory[i]);

        // The FFT class provides no FFT normalization convention between fwd and rev transf. Need to apply the correct one.
        double T = (double)v_tilde->size();
        double sqrT = sqrt(T);
        for(size_t n=0; n<v_tilde->size(); n++)
        {
            (*v_tilde)[n].re_/= sqrT;
            (*v_tilde)[n].im_/= sqrT;
        }

        // Calculate normalization
        double C = 0.0;
        for(size_t t=0; t<atomVelTrajectory[i].size(); t++)
        {
            double v = atomVelTrajectory[i][t].re_;
            C+= v*v;
        }

        // Compute the normalized dot products between the fft and its complex conjugates,
        // i.e. $\mathbf{\tilde{v}}^{*}_n \cdot \mathbf{\tilde{v}}_n / (C \Delta \omega)$
        double T2 = T*T;
        double dw = 2.0 * M_PI * (T - 1.0) / (timeStep * T2);
        double dwC =  C * dw;
        for(int n=0; n<v_tilde->size(); n++)
        {
            avg_v_tilde2[n]+= ( (*v_tilde)[n].modulus2() / dwC );
        }

        pb.updateProgress(i, (int)atomIndicesToInclude.size());
    }

    pb.endProgress();


    // Average over atoms
    for(int n=0; n<avg_v_tilde2.size(); n++)
    {
        avg_v_tilde2[n]/= double(numAtoms);
    }


    // Print the resulting VDOS spectrum
    int T = (int)avg_v_tilde2.size();
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s%-15s\r\n", "omega[rad/fs]", "vdos");
    for(int n=0; n<T; n++)
    {
        double omega = 2.0*M_PI*double(n) / (double(T)*timeStep);
        fprintf(stdOut_, "\t% -15.8f% -15.8f\r\n", omega, avg_v_tilde2[n]);
    }

    return lastError_;
}
