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

#include "CmdFromXTC.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../../Utilities/XTCFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdFromXTC::getCmd()
{
    return "fromxtc";
}

std::vector<std::string> CCmdFromXTC::getCmdLineKeywords()
{
    return { "fromxtc" };
}

std::vector<std::string> CCmdFromXTC::getCmdHelpLines()
{
    return {
                "fromxtc <XTC filename> <DCD filename> <Num. time steps>"
           };
}

std::string CCmdFromXTC::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tConverts an XTC file, <XTC filename>, to a DCD file with name <DCD filename>.\r\n";
    text+= "\tThe number of time steps is not available in the XTC header. Hence, this must\r\n";
    text+= "\tbe given in the <Num. time steps> argument (i.e., the total number of time\r\n";
    text+= "\tsteps in the simulation).";

    return text;
}

std::string CCmdFromXTC::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    #if GROMACS_LIB_INSTALLED == 1
        #ifdef GROMACS_XTC_HEADERS_INSTALLED
            size_t arg = 0;
            CProgressBar pb;

            std::string xtcFileName = CASCIIUtility::getArg(arguments, arg++);
            std::string dcdFileName = CASCIIUtility::getArg(arguments, arg++);
            int numTimeSteps = atoi(CASCIIUtility::getArg(arguments, arg++).data());

            if(numTimeSteps <= 0)
            {
                lastError_ = "There must be a positive number of timesteps!";
                return lastError_;
            }

            // Open XTC file first to calculate stride, etc.
            t_fileio* file = open_xtc(xtcFileName.c_str(), "r");
            if(!file)
            {
                lastError_ = std::string("Could not open file: ") + xtcFileName + std::string("!");
                return lastError_;
            }

            // Set up variables to receive a single frame
            int numAtoms = 0;
            int64_t step = 0;
            real time = 0;
            matrix box = {{ 1.0f, 0.0f, 0.0f }, { 0.0f, 1.0f, 0.0f }, { 0.0f, 0.0f, 1.0f }};
            rvec* coordinates;
            real precision = 0.0f;
            gmx_bool xtcOK = false;

            // Read first and second entry. Gather info and then close
            read_first_xtc(file, &numAtoms, &step, &time, box, &coordinates, &precision, &xtcOK);
            if(!xtcOK)
            {
                lastError_ = std::string("Could not read first XTC entry in ") + xtcFileName + std::string("!");
                close_xtc(file);
                return lastError_;
            }
            int step1 = (int)step;
            float time1 = (float)time;

            read_next_xtc(file, numAtoms, &step, &time, box, coordinates, &precision, &xtcOK);
            if(!xtcOK)
            {
                lastError_ = std::string("Could not read second XTC entry in ") + xtcFileName + std::string("!");
                close_xtc(file);
                return lastError_;
            }
            int outputStride = step - step1;
            float timeStep = (outputStride != 0) ? (time - time1) / float(outputStride) : (time - time1);

            close_xtc(file);

            // Create a lambda to create a DCD frame from XTC coordinates
            std::function<void(const std::string&, const int&, const matrix&, const rvec*)> createDCDFrameFromXTCCoord =
                    [](const std::string& dcdFileName, const int& numAtoms, const matrix& box, const rvec* coordinates)
            {
                std::function<std::tuple<double, double, double>(const int&)> getAtomPos = [coordinates](const int& atomIndex)
                {
                    double X = (double)coordinates[atomIndex][0] * 10.0f; // Convert from nm to AA
                    double Y = (double)coordinates[atomIndex][1] * 10.0f; // Convert from nm to AA
                    double Z = (double)coordinates[atomIndex][2] * 10.0f; // Convert from nm to AA

                    return std::tuple<double, double, double>(X, Y, Z);
                };

                CDCDFile::appendToDCDFile(dcdFileName, numAtoms, { box[0][0]*10.0f, box[1][1]*10.0f, box[2][2]*10.0f }, getAtomPos);
            };

            // Open XTC file again and read firs entry for conversion
            file = open_xtc(xtcFileName.c_str(), "r");
            read_first_xtc(file, &numAtoms, &step, &time, box, &coordinates, &precision, &xtcOK);
            CDCDFile::createDCDFileIfNotExists(dcdFileName, numTimeSteps, outputStride, timeStep, numAtoms);
            createDCDFrameFromXTCCoord(dcdFileName, numAtoms, box, coordinates);

            // Read remaining frames from XTC and convert to DCD
            int64_t prevStep = step;
            pb.beginProgress("Converting XTC file contents to DCD");
            while(1)
            {
                read_next_xtc(file, numAtoms, &step, &time, box, coordinates, &precision, &xtcOK);
                if(xtcOK && (step != prevStep))
                {
                    createDCDFrameFromXTCCoord(dcdFileName, numAtoms, box, coordinates);
                    pb.updateProgress(step % 1000, 1000);
                    prevStep = step;
                }
                else
                {
                    break;
                }
            }

            pb.endProgress();
            close_xtc(file);

        #else
            lastError_ = "Error, missing headers: MolTwister was not built with the Gromacs XTC library, likely because it was not available during the MolTwister build. Could not load XTC file!";
            return lastError_;
        #endif
    #else
        lastError_ = "Error, missing library. MolTwister was not built with the Gromacs XTC library, likely because it was not available during the MolTwister build. Could not load XTC file!";
        return lastError_;
    #endif

    return lastError_;
}
