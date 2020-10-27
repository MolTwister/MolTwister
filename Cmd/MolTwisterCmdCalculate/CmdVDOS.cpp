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
                "vdos <DCD filename> <frame from> <frame to> <time step (fs)> name <atom IDs (comma sep., no space)>",
                "vdos <DCD filename> <frame from> <frame to> <time step (fs)> sel"
           };
}

std::string CCmdVDOS::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the vibrational density of states (VDOS) for the DCD file, <DCD filename>,\r\n";
    text+= "\tbetween the start frame, <frame from>, to the frame <frame to>. The time step, <time\r\n";
    text+= "\tstep>, given in units of femtoseconds (fs), is used to obtain numerically calculated\r\n";
    text+= "\tvelocities that are needed within the calculations. The VDOS is only calculated for a\r\n";
    text+= "\tgiven selection of atoms. Either, based on the atom 'name', where a list, <atom IDs>,\r\n";
    text+= "\tis supplied (e.g., H, O, C7), or based on the visual selection of atoms, achieved through\r\n";
    text+= "\tthe 'sel' keyword.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. omega vdos\r\n";
    text+= "\t2. <angular frequency> <VDOS value>\r\n";
    text+= "\t3. <angular frequency> <VDOS value>\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\tN+1. <angular frequency> <VDOS value>\r\n";
    text+= "\twhere N is the length of the VDOS function.";

    return text;
}

std::string CCmdVDOS::execute(std::vector<std::string> arguments)
{
    /////////////////////////////////////////////////////////////////////////////////
    // Calculate the vibrational density of states (VDOS) by Fourier transformation
    // of the VACF (which in turn can be found through FFT of the VDOS).
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


    // Pre-read all selected atom trajectories
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
            atomTrajectory[n][i].resize(fftLen);
        }
    }

    pb. beginProgress("Reading atomic trajectories");
    for(int t=frameFrom; t<frameTo; t++)
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
            atomTrajectory[0][i][t].re_ = v.x_;
            atomTrajectory[1][i][t].re_ = v.y_;
            atomTrajectory[2][i][t].re_ = v.z_;
        }

        pb.updateProgress(t - frameFrom, fftLen);
    }
    pb.endProgress();


    // Loop through each particle, then FFT each of the particle trajectories
    std::vector<double> avg_v_tilde2;
    avg_v_tilde2.resize(fftLen, 0.0);
    pb.beginProgress("Calculating vibrational density of states");
    for(int i=0; i<numAtoms; i++)
    {
        // Fourier transform trajectory of atom i in x,y and z directions
        CFFT1D fft;
        std::shared_ptr<std::vector<CFFT1D::CCplx>> v_tilde[3];
        v_tilde[0] = fft.fft1D(atomTrajectory[0][i]);
        v_tilde[1] = fft.fft1D(atomTrajectory[1][i]);
        v_tilde[2] = fft.fft1D(atomTrajectory[2][i]);

        // Compute the dot products between the fft and its complex conjugates,
        // i.e. $\mathbf{\tilde{v}}^{*}_n \cdot \mathbf{\tilde{v}}_n$
        for(int n=0; n<v_tilde[0]->size(); n++)
        {
            double vx2 = (*v_tilde[0])[n].re_*(*v_tilde[0])[n].re_ + (*v_tilde[0])[n].im_*(*v_tilde[0])[n].im_;
            double vy2 = (*v_tilde[1])[n].re_*(*v_tilde[1])[n].re_ + (*v_tilde[1])[n].im_*(*v_tilde[1])[n].im_;
            double vz2 = (*v_tilde[2])[n].re_*(*v_tilde[2])[n].re_ + (*v_tilde[2])[n].im_*(*v_tilde[2])[n].im_;
            avg_v_tilde2[n]+= (vx2 + vy2 + vz2);
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
    int M = (int)avg_v_tilde2.size();
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s%-15s\r\n", "omega", "vdos");
    for(int n=0; n<M; n++)
    {
        double omega = 2.0*M_PI*double(n) / (double(M)*timeStep);
        fprintf(stdOut_, "\t% -15.8f% -15.8f\r\n", omega, avg_v_tilde2[n]);
    }

    return lastError_;
}
