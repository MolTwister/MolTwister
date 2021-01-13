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

#include "CmdVACF.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdVACF::getCmd()
{
    return "vacf";
}

std::vector<std::string> CCmdVACF::getCmdLineKeywords()
{
    return { "vacf", "name", "sel" };
}

std::vector<std::string> CCmdVACF::getCmdHelpLines()
{
    return {
                "vacf <DCD filenmae> <frame from> <frame to> <time step (fs)> <VACF length> name <atom IDs (comma sep., no space)>",
                "vacf <DCD filenmae> <frame from> <frame to> <time step (fs)> <VACF length> sel"
           };
}

std::string CCmdVACF::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the velocity auto correlation function (VACF) for the DCD file, <DCD filename>,\r\n";
    text+= "\tbetween the start frame, <frame from>, to the frame <frame to>. The time step, <time step>,\r\n";
    text+= "\tgiven in units of femtoseconds (fs), is used to obtain numerically calculated velocities\r\n";
    text+= "\tfrom pairs of atomic position originating from two adjacent time frames. The <VACF length>\r\n";
    text+= "\tparameter is the length of the VACF, in number of time steps, and must be smaller or equal\r\n";
    text+= "\tto D = <frame to> - <frame from>. If <VACF length> < D, then a VACF is calculated for each\r\n";
    text+= "\tstarting time t0 from t0=<frame from> to t0=<frame to> (where indices outside the number\r\n";
    text+= "\tof frames are ignored). All collected VACF are then averaged. The VACF is only calculated\r\n";
    text+= "\tfor a given selection of atoms. Either, based on the atom 'name', where a list, <atom IDs>,\r\n";
    text+= "\tis supplied (e.g., H, O, C7), or based on the visual selection of atoms, achieved through\r\n";
    text+= "\tthe 'sel' keyword.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. t vacf norm(vacf)\r\n";
    text+= "\t2. <VACF time point in fs> <VACF value> <VACF value, normalized to the first value>\r\n";
    text+= "\t3. <VACF time point in fs> <VACF value> <VACF value, normalized to the first value>\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\tN+1. <VACF time point in fs> <VACF value> <VACF value, normalized to the first value>\r\n";
    text+= "\twhere N is the number of points selected in the VACF (<VACF length>).";

    return text;
}

std::string CCmdVACF::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<double> vacf;
    std::vector<int> vacfAvgCount;
    std::vector<int> aAtomIndicesToInclude;
    CDCDFile dcdFile;
    std::string text;
    double timeStep;
    int numFrames, frameFrom, frameTo, vacfLen;


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
    frameTo = atoi(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    timeStep = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vacfLen = atoi(text.data());
    if(timeStep < 0.0)
    {
        lastError_ = "Error: timestep cannot be zero!";
        return lastError_;
    }
    if(vacfLen < 1)
    {
        lastError_ = "Error: vacf length must be greater than zero!";
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
                    aAtomIndicesToInclude.emplace_back(i);
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
                aAtomIndicesToInclude.emplace_back(i);
            }
        }
    }
    else
    {
        lastError_ = "Error: expected 'name' or 'sel'!";
        return lastError_;
    }


    // Loop through different starting times, i0 (i.e. t0)
    vacf.resize(vacfLen, 0.0);
    vacfAvgCount.resize(vacfLen, 0);

    pb.beginProgress("Calculating vacf");
    for(int i0=0; i0<(frameTo-frameFrom); i0++)
    {
        // Loop through the frames starting from i0 (i.e. t0)
        bool firstPosRead = false;
        bool initialVelocityCalculated = false;
        std::vector<C3DVector> r1;
        std::vector<C3DVector> r2;
        std::vector<C3DVector> v0, vi;
        numFrames = dcdFile.getNumRecords();
        for(int i=(frameFrom+i0); i<(frameFrom+i0+vacfLen); i++)
        {
            if((i < 0) || (i >= numFrames)) continue;

            dcdFile.gotoRecord(i);

            // If we have r1 we can proceed to calculate r2 and velocities
            // from r1,r2, else we must first store r1.
            if(firstPosRead)
            {
                // Obtain r2 from all selected atoms
                r2.clear();
                for(int j=0; j<aAtomIndicesToInclude.size(); j++)
                {
                    r2.emplace_back(dcdFile.getCoordinate(aAtomIndicesToInclude[j]));
                }

                // If we have not calculated v0, we must do that now (and assign vi=v0), else
                // we can proceed to calculate vi, which can be used for <v(t0)*v(t0+t)>
                if(!initialVelocityCalculated)
                {
                    v0.clear();
                    for(int j=0; j<r2.size(); j++)
                    {
                        C3DVector vel;

                        if(j >= r1.size()) continue;
                        vel.x_ = (r2[j].x_ - r1[j].x_) / timeStep; // [AA/fs]
                        vel.y_ = (r2[j].y_ - r1[j].y_) / timeStep; // [AA/fs]
                        vel.z_ = (r2[j].z_ - r1[j].z_) / timeStep; // [AA/fs]

                        v0.emplace_back(vel);
                    }

                    vi = v0;
                    initialVelocityCalculated = true;
                }
                else
                {
                    vi.clear();
                    for(int j=0; j<r2.size(); j++)
                    {
                        C3DVector vel;

                        if(j >= r1.size()) continue;
                        vel.x_ = (r2[j].x_ - r1[j].x_) / timeStep; // [AA/fs]
                        vel.y_ = (r2[j].y_ - r1[j].y_) / timeStep; // [AA/fs]
                        vel.z_ = (r2[j].z_ - r1[j].z_) / timeStep; // [AA/fs]

                        vi.emplace_back(vel);
                    }
                }

                // Calculate <v(t0)*v(t0+t)> now that we have both v(t0) and v(t0+t)
                double VACF = 0.0;
                for(int j=0; j<v0.size(); j++)
                {
                    if(j >= vi.size()) continue;
                    VACF+= v0[j]*vi[j];
                }
                VACF/= double(v0.size());

                // Store the result
                vacf[i-(frameFrom+i0)]+= VACF;
                vacfAvgCount[i-(frameFrom+i0)]++;

                // Make r2 the next r1
                r1.clear();
                r1 = r2;
            }
            else
            {
                r1.clear();
                for(int j=0; j<aAtomIndicesToInclude.size(); j++)
                {
                    r1.emplace_back(dcdFile.getCoordinate(aAtomIndicesToInclude[j]));
                }
                firstPosRead = true;
            }

            pb.updateProgress(i0, frameTo-frameFrom);
        }
    }
    pb.endProgress();


    // Calculate the first value of vacf used for normalization
    double nvacf;
    if(vacf.size() > 1) nvacf = vacf[1] / double(vacfAvgCount[1]);
    else                nvacf = 1.0;


    // Print the results
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s%-15s%-15s\r\n", "t", "vacf", "norm(vacf)");
    for(int i=1; i<vacf.size(); i++)
    {
        double VACF = vacf[i] / double(vacfAvgCount[i]);
        fprintf(stdOut_, "\t% -15.8f% -15.8f% -15.8f\r\n", double(i-1)*timeStep, VACF, VACF / nvacf);
    }

    return lastError_;
}
