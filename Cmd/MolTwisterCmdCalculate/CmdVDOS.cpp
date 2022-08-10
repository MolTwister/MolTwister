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

#include "CmdVDOS.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdVDOS::getCmd()
{
    return "vdos";
}

std::vector<std::string> CCmdVDOS::getCmdLineKeywords()
{
    return { "vdos", "name", "sel", "com" };
}

std::vector<std::string> CCmdVDOS::getCmdHelpLines()
{
    return {
                "vdos <DCD filename> <frame from> <num. bits> <time step (fs)> [com] name <atom IDs (comma sep., no space)>",
                "vdos <DCD filename> <frame from> <num. bits> <time step (fs)> [com] sel"
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
    text+= "\tthe 'sel' keyword. If the 'com' keyword is selected, the velocity of the center of mass of\r\n";
    text+= "\tthe selected atoms is calculated, otherwise the average velocity of the selected atoms is used.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. omega[rad/fs] vdos_o freq[x10^15Haz] vdos_f lambda[nm] vdos_l k[cm^-1] vdos_k\r\n";
    text+= "\t2. <angular frequency> <VDOS value> <frequency> <VDOS value> <wavelength> <VDOS value> <wavenumber> <VDOS value>\r\n";
    text+= "\t3. <angular frequency> <VDOS value> <frequency> <VDOS value> <wavelength> <VDOS value> <wavenumber> <VDOS value>\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\tN+1. <angular frequency> <VDOS value> <frequency> <VDOS value> <wavelength> <VDOS value> <wavenumber> <VDOS value>\r\n";
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


    // Check if we are going to use COM or average (default is average)
    bool useCOMofAtoms = false;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "com")
    {
        useCOMofAtoms = true;
    }
    else arg--;


    // Pre-read all selected atom trajectories (that need to be fftLen+1, while velocity traj. becomes fftLen long)
    int fftLen = frameTo - frameFrom;
    int numAtoms = (int)atomIndicesToInclude.size();

    std::vector<std::vector<std::vector<CFFT1D::CCplx>>> atomTrajectory;
    if(useCOMofAtoms)
    {
        lastError_ = fillAtomTrajectoryCOMOfAtoms(atomIndicesToInclude, frameFrom, frameTo, fftLen, atomTrajectory, pb, dcdFile);
    }
    else
    {
        lastError_ = fillAtomTrajectoryAvgAtoms(atomIndicesToInclude, frameFrom, frameTo, fftLen, atomTrajectory, pb, dcdFile);
    }
    if(!lastError_.empty()) return lastError_;


    // Convert positional trajectories to velocity trajectories
    pb.beginProgress("Reading atomic trajectories");
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

    std::vector<std::vector<CFFT1D::CCplx>> atomVelTrajectory[3] = { std::vector<std::vector<CFFT1D::CCplx>>(trajAtomCount),
                                                                     std::vector<std::vector<CFFT1D::CCplx>>(trajAtomCount),
                                                                     std::vector<std::vector<CFFT1D::CCplx>>(trajAtomCount) };
    for(size_t i=0; i<trajAtomCount; i++)
    {
        for(int k=0; k<3; k++) atomVelTrajectory[k][i].resize(fftLen);
    }

    for(size_t t=1; t<trajSize; t++)
    {
        for(size_t i=0; i<trajAtomCount; i++)
        {
            for(int k=0; k<3; k++)
            {
                double d1 = atomTrajectory[k][i][t-1].re_;
                double d2 = atomTrajectory[k][i][t].re_;
                double v = (d2 - d1) / timeStep; // AA / fs
                atomVelTrajectory[k][i][t-1].re_ = v;
            }
        }

        pb.updateProgress((int)t, (int)trajSize);
    }

    pb.endProgress();


    // Loop through each particle, then FFT each of the particle velocity trajectories
    std::vector<double> avg_v_tilde2[3];
    for(int k=0; k<3; k++) avg_v_tilde2[k].resize(fftLen, 0.0);
    pb.beginProgress("Calculating vibrational density of states");
    for(int i=0; i<(int)trajAtomCount; i++)
    {
        // Calculate normalization
        double C = 0.0;
        for(int k=0; k<3; k++)
        {
            for(size_t t=0; t<atomVelTrajectory[k][i].size(); t++)
            {
                double v = atomVelTrajectory[k][i][t].re_;
                C+= v*v;
            }
        }

        for(int k=0; k<3; k++)
        {
            // Fourier transform velocity trajectory of atom i
            CFFT1D fft;
            std::shared_ptr<std::vector<CFFT1D::CCplx>> v_tilde;
            v_tilde = fft.fft1D(atomVelTrajectory[k][i]);

            // The FFT class provides no FFT normalization convention between fwd and rev transf. Need to apply the correct one.
            double T = (double)v_tilde->size();
            double sqrT = sqrt(T);
            for(size_t n=0; n<v_tilde->size(); n++)
            {
                (*v_tilde)[n].re_/= sqrT;
                (*v_tilde)[n].im_/= sqrT;
            }

            // Compute the normalized dot products between the fft and its complex conjugates,
            // i.e. $\mathbf{\tilde{v}}^{*}_n \cdot \mathbf{\tilde{v}}_n / (C \Delta \omega)$
            double T2 = T * T;
            double dw = 2.0 * M_PI * (T - 1.0) / (timeStep * T2);
            double dwC =  C * dw;
            for(int n=0; n<v_tilde->size(); n++)
            {
                avg_v_tilde2[k][n]+= ( (*v_tilde)[n].modulus2() / dwC );
            }
        }

        pb.updateProgress(i, (int)atomIndicesToInclude.size());
    }

    pb.endProgress();


    // Average over atoms
    for(int k=0; k<3; k++)
    {
        for(int n=0; n<avg_v_tilde2[k].size(); n++)
        {
            avg_v_tilde2[k][n]/= double(numAtoms);
        }
    }


    // Sum v_x^2 + v_y^2 + v_z^2 to form dot product and store result in v_x^2 vector
    for(int k=1; k<3; k++)
    {
        for(int n=0; n<avg_v_tilde2[k].size(); n++)
        {
            avg_v_tilde2[0][n]+= avg_v_tilde2[k][n];
        }
    }


    // Print the resulting VDOS spectrum
    int T = (int)avg_v_tilde2[0].size();
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\r\n", "omega[rad/fs]", "vdos_o", "f[x10^15Hz]", "vdos_f", "lambda[nm]", "vdos_l", "k[cm^-1]", "vdos_k");
    for(int n=0; n<T; n++)
    {
        double freq = double(n) / (double(T)*timeStep*1.0E-15); // [Hz]
        double omega = 2.0 * M_PI * freq; // [rad/s]
        double lambda = (3.0E8 / freq); // [m]
        double k = omega / (2.0*M_PI * 3.0E8); // [m^-1]
        fprintf(stdOut_, "\t% -15.8f% -15.8f% -15.8f% -15.8f% -15.8f% -15.8f% -15.8f% -15.8f\r\n",
                omega*1.0E-15, // omega
                avg_v_tilde2[0][n], // vdos_o
                freq*1.0E-15, // f
                2.0*M_PI*avg_v_tilde2[0][n], // vdos_f
                lambda*1.0E9, // lambda
                avg_v_tilde2[0][n] * (2.0*M_PI*3.0E8 / (lambda*lambda)) * 1.0E-24, // vdos_o ( 1/(m*s) = 1/(nm*fs*10^24) )
                k*1.0E-2, // k
                avg_v_tilde2[0][n] * 2.0*M_PI * 0.00003); // vdos_k ( c = 3E8 m/s = 3E8 * 10^2 * 10^-15 cm/fs = 0.00003 cm/fs )
    }

    return lastError_;
}

std::string CCmdVDOS::fillAtomTrajectoryAvgAtoms(const std::vector<int>& atomIndicesToInclude, int frameFrom, int frameTo, int fftLen, std::vector<std::vector<std::vector<CFFT1D::CCplx>>>& atomTrajectory, CProgressBar& pb, CDCDFile &dcdFile) const
{
    int numAtoms = (int)atomIndicesToInclude.size();
    int numFrames = dcdFile.getNumRecords();

    atomTrajectory.resize(3);
    for(int n=0; n<3; n++)
    {
        atomTrajectory[n].resize(numAtoms);
        for(int i=0; i<numAtoms; i++)
        {
            atomTrajectory[n][i].resize(fftLen+1);
        }
    }

    pb.beginProgress("Reading atomic trajectories");
    for(int t=frameFrom; t<(frameTo+1); t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            return std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
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

    return "";
}

std::string CCmdVDOS::fillAtomTrajectoryCOMOfAtoms(const std::vector<int>& atomIndicesToInclude, int frameFrom, int frameTo, int fftLen, std::vector<std::vector<std::vector<CFFT1D::CCplx>>>& atomTrajectory, CProgressBar& pb, CDCDFile& dcdFile) const
{
    int numAtoms = (int)atomIndicesToInclude.size();
    int numFrames = dcdFile.getNumRecords();

    atomTrajectory.resize(3);
    for(int n=0; n<3; n++)
    {
        atomTrajectory[n].resize(1);
        for(int i=0; i<numAtoms; i++)
        {
            atomTrajectory[n][i].resize(fftLen+1);
        }
    }

    pb.beginProgress("Reading atomic trajectories");
    for(int t=frameFrom; t<(frameTo+1); t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            return std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
        }

        C3DVector v;
        dcdFile.gotoRecord(t);
        double M = 0.0;
        for(int i=0; i<numAtoms; i++)
        {
            int atomIndex = atomIndicesToInclude[i];

            double m = state_->atoms_[atomIndex]->m_;
            M+= m;

            v = dcdFile.getCoordinate(atomIndex);
            atomTrajectory[0][0][t-frameFrom].re_+= (m * v.x_);
            atomTrajectory[1][0][t-frameFrom].re_+= (m * v.y_);
            atomTrajectory[2][0][t-frameFrom].re_+= (m * v.z_);
        }

        atomTrajectory[0][0][t-frameFrom].re_/= M;
        atomTrajectory[1][0][t-frameFrom].re_/= M;
        atomTrajectory[2][0][t-frameFrom].re_/= M;

        pb.updateProgress(t-frameFrom, fftLen+1);
    }
    pb.endProgress();

    return "";
}
