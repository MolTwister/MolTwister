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

#include "CmdDistProbabilityCOM.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdDistProbabilityCOM::getCmd()
{
    return "distprobabilitycom";
}

std::vector<std::string> CCmdDistProbabilityCOM::getCmdLineKeywords()
{
    return { "distprobabilitycom", "indices", "atomtypes", "ignorepbc" };
}

std::vector<std::string> CCmdDistProbabilityCOM::getCmdHelpLines()
{
    return {
                "distprobabilitycom <DCD filename> <frame from> <frame to> <num dist. bins> <start dist. range> <end dist. range> <num COM bins> <COM direction> <start COM range> <end COM range> <atom IDs from> <atom IDs to> indices <M> <COM base index 1> ... <COM base index M> [ignorepbc]",
                "distprobabilitycom <DCD filename> <frame from> <frame to> <num dist. bins> <start dist. range> <end dist. range> <num COM bins> <COM direction> <start COM range> <end COM range> <atom IDs from> <atom IDs to> atomtypes <M> <COM base ID 1> ... <COM base ID M> [ignorepbc]"
           };
}

std::string CCmdDistProbabilityCOM::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates a two dimensional surface plot with the molecular center of mass (COM) along\r\n";
    text+= "\tone axis. The other axis is the length between atoms of the set defined by <atom IDs from>\r\n";
    text+= "\t(comma separated, no space) and the set defined by <atom IDs to>. Each of the axis are\r\n";
    text+= "\tbinned according to <num COM bins> and <num dist. bins>, respectively. The range of the\r\n";
    text+= "\ttwo axis are defined by [<start COM range>, <end COM range>] and by [<start dist. range>,\r\n";
    text+= "\t<end dist. range>], respectively. The <COM direction> can be either x, y or z and defines\r\n";
    text+= "\twhich of the components of the COM vector are to be binned. The set of atoms that is used\r\n";
    text+= "\tas basis for the COM calculations are defined by specifying the number of atoms, <M>,\r\n";
    text+= "\tfollowed by the atom indices <COM base index i>, if 'indices' is specified, or by the\r\n";
    text+= "\tatom IDs (e.g., H, O, C7) <COM base ID i>, if 'atomtypes' is specified. By default the\r\n";
    text+= "\tperiodic boundary conditions (PBC) are taken into account (i.e., distances are measured\r\n";
    text+= "\tacross PBC boundaries. However, PBCs can be ignored by using the 'ignorepbc' keyword.\r\n";
    text+= "\t\r\n";
    text+= "\tTrajectories used as input is specified by <DCD filename>, where all frames from <frame from>\r\n";
    text+= "\tto <frame to> are used in the averaging (i.e., contributing to counts in the 2D grid of\r\n";
    text+= "\tspecified bins).\r\n";
    text+= "\t\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Num distance bins: <num dist. bins>\r\n";
    text+= "\t2. Start dist. range: <start dist. range>\r\n";
    text+= "\t3. End dist. range: <end dist. range>\r\n";
    text+= "\t4. Num COM bins: <num COM bins>\r\n";
    text+= "\t5. Start COM range: <start COM range>\r\n";
    text+= "\t6. End COM range: <end COM range>\r\n";
    text+= "\t7. COM range direction (0=x,1=y,2=z): <COM direction>\r\n";
    text+= "\t8. From frame: <frame from>\r\n";
    text+= "\t9. To frame: <frame to>\r\n";
    text+= "\t10. Dist from: <atom IDs from>\r\n";
    text+= "\t11. Dist to: <atom IDs to>\r\n";
    text+= "\t12. COMPos Dist Tot.Count Normalized\r\n";
    text+= "\t13. <COM positio binn> <distance between atoms bin> <total count in bin> <normalized count, max set to unity>\r\n";
    text+= "\t14. <COM positio binn> <distance between atoms bin> <total count in bin> <normalized count, max set to unity>\r\n";
    text+= "\t                .\r\n";
    text+= "\t                .\r\n";
    text+= "\t                .\r\n";
    text+= "\tNa*Nb+12. <COM positio binn> <distance between atoms bin> <total count in bin> <normalized count, max set to unity>\r\n";
    text+= "\twhere Na and Nb are the number of COM bins and distance bins, respectively.";

    return text;
}

std::string CCmdDistProbabilityCOM::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToIncludeDistFrom;
    std::vector<int> atomIndicesToIncludeDistTo;
    std::vector<int> atomIndicesToIncludeCOM;
    CDCDFile dcdFile;
    std::string text;
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
    frameTo = atoi(text.data());

    if((frameTo - frameFrom) <= 0)
    {
        lastError_ = "Error: no frames were selected!";
        return lastError_;
    }


    // Retrieve the number of bins to divide dist. range into
    text = CASCIIUtility::getArg(arguments, arg++);
    int numBins = atoi(text.data());


    // Retrieve the dist. range
    text = CASCIIUtility::getArg(arguments, arg++);
    double startRange = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double endRange = atof(text.data());

    if(fabs(startRange - endRange) < double(FLT_MIN))
    {
        lastError_ = "Error: cannot have a start range equal to end range!";
        return lastError_;
    }


    // Retrieve the number of bins to divide com range into
    text = CASCIIUtility::getArg(arguments, arg++);
    int numCOMBins = atoi(text.data());


    // Retrieve the COM direction to bin
    int comDir = -1;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "x")         comDir = 0;
    else if(text == "y")    comDir = 1;
    else if(text == "z")    comDir = 2;
    else
    {
        lastError_ = "Error: expected 'x', 'y' or 'z'!";
        return lastError_;
    }


    // Retrieve the COM range in previously chosen direction
    text = CASCIIUtility::getArg(arguments, arg++);
    double startRangeCOM = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double endRangeCOM = atof(text.data());

    if(fabs(startRangeCOM - endRangeCOM) < double(FLT_MIN))
    {
        lastError_ = "Error: cannot have a COM start range equal to end range!";
        return lastError_;
    }


    // Retrieve the atoms to calculate distance from
    std::string distFromString;
    distFromString = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(distFromString);
    atomsToInclude = CASCIIUtility::getWords(distFromString, ",");

    for(int i=0; i<atomsToInclude.size(); i++)
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            std::string ID = state_->atoms_[j]->getID();
            if(ID == atomsToInclude[i])
            {
                atomIndicesToIncludeDistFrom.emplace_back(j);
            }
        }
    }


    // Retrieve the atoms to calculate distance to
    std::string distToString;
    distToString = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(distToString);
    atomsToInclude = CASCIIUtility::getWords(distToString, ",");

    for(int i=0; i<atomsToInclude.size(); i++)
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            std::string ID = state_->atoms_[j]->getID();
            if(ID == atomsToInclude[i])
            {
                atomIndicesToIncludeDistTo.emplace_back(j);
            }
        }
    }


    // Retrieve atom types for which to calculate COM
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "indices")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        int numIndices = atoi(text.data());

        for(int i=0; i<numIndices; i++)
        {
            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndicesToIncludeCOM.emplace_back(atoi(text.data()));
        }
    }
    else if(text == "atomtypes")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        int numAtomTypes = atoi(text.data());

        for(int i=0; i<numAtomTypes; i++)
        {
            text = CASCIIUtility::getArg(arguments, arg++);

            for(int j=0; j<state_->atoms_.size(); j++)
            {
                std::string ID = state_->atoms_[j]->getID();
                if(ID == text)
                {
                    atomIndicesToIncludeCOM.emplace_back(j);
                }
            }
        }
    }


    // Retrieve PBC
    if(state_->currentFrame_ < 0)
    {
        lastError_ = "Error: invalid frame!";
        return lastError_;
    }
    if(!state_->view3D_)
    {
        lastError_ = "Error: could not communicate with 3D view!";
        return lastError_;
    }
    C3DRect* pbcPtr = nullptr;
    C3DRect pbc = state_->view3D_->calcPBC(state_->currentFrame_);
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text != "ignorepbc")
    {
        pbcPtr = &pbc;
        arg--;
    }


    // Loop through all distances to measure and bin them
    int numProcessedFrames = 0;
    int numFrames = dcdFile.getNumRecords();
    std::vector<std::vector<int>> comDistBins;
    std::vector<std::vector<double>> comDistBinsNormalized;
    comDistBins.resize(numCOMBins);
    comDistBinsNormalized.resize(numCOMBins);
    for(int i=0; i<comDistBins.size(); i++)
    {
        comDistBins[i].resize(numBins, 0);
        comDistBinsNormalized[i].resize(numBins, 0.0);
    }

    pb.beginProgress("Calculating distance-COM distributions");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);
        CDCDTools dcdTools(state_, stdOut_);
        C3DVector vecCOM = dcdTools.getCenterOfMass(atomIndicesToIncludeCOM, &dcdFile);
        double com = 0.0;
        if(comDir == 0) com = vecCOM.x_;
        if(comDir == 1) com = vecCOM.y_;
        if(comDir == 2) com = vecCOM.z_;

        for(int i=0; i<atomIndicesToIncludeDistFrom.size(); i++)
        {
            for(int j=0; j<atomIndicesToIncludeDistTo.size(); j++)
            {
                double dDist = dcdTools.measureDist(dcdFile, atomIndicesToIncludeDistFrom[i], atomIndicesToIncludeDistTo[j], pbcPtr);
                int iBin = int((dDist - startRange) / (endRange - startRange) * double(numBins));
                int iBinCOM = int((com - startRangeCOM) / (endRangeCOM - startRangeCOM) * double(numCOMBins));

                if((iBin >= 0) && (iBin < numBins) && (iBinCOM >= 0) && (iBinCOM < numCOMBins)) comDistBins[iBinCOM][iBin]++;
            }
        }

        numProcessedFrames++;
        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Normalize the distance distribution in each bin, where
    // maximum is set to unity
    for(int i=0; i<comDistBins.size(); i++)
    {
        int max = 0;

        // Find maximum
        for(int j=0; j<comDistBins[i].size(); j++)
        {
            if(comDistBins[i][j] > max) max = comDistBins[i][j];
        }

        // Create normalized value
        for(int j=0; j<comDistBins[i].size(); j++)
        {
            if(max > 0)    comDistBinsNormalized[i][j] = double(comDistBins[i][j]) / double(max);
            else           comDistBinsNormalized[i][j] = 0.0;
        }
    }


    // Present the results
    printf("\r\n");
    fprintf(stdOut_, "\tNum distance bins: %i\r\n", numBins);
    fprintf(stdOut_, "\tStart dist. range: %g\r\n", startRange);
    fprintf(stdOut_, "\tEnd dist. range: %g\r\n", endRange);
    fprintf(stdOut_, "\tNum COM bins: %i\r\n", numCOMBins);
    fprintf(stdOut_, "\tStart COM range: %g\r\n", startRangeCOM);
    fprintf(stdOut_, "\tEnd COM range: %g\r\n", endRangeCOM);
    fprintf(stdOut_, "\tCOM range direction (0=x,1=y,2=z): %i\r\n", comDir);
    fprintf(stdOut_, "\tFrom frame: %i\r\n", frameFrom);
    fprintf(stdOut_, "\tTo frame: %i\r\n", frameTo);
    fprintf(stdOut_, "\tDist from: %s\r\n\r\n", distFromString.data());
    fprintf(stdOut_, "\tDist to: %s\r\n\r\n", distToString.data());
    fprintf(stdOut_, "\t%-15s%-15s%-15s%-15s\r\n", "COMPos", "Dist", "Tot.Count", "Normalized");
    for(int i=0; i<numCOMBins; i++)
    {
        for(int j=0; j<numBins; j++)
        {
            fprintf(stdOut_, "\t% -15.4g% -15.4g%-15i%-15.3f\r\n", 0.5 * (endRangeCOM - startRangeCOM) / double(numCOMBins) * (2.0*double(i) + 1.0) + startRangeCOM,
                                                                   0.5 * (endRange - startRange) / double(numBins) * (2.0*double(j) + 1.0) + startRange,
                                                                   comDistBins[i][j],
                                                                   comDistBinsNormalized[i][j]);
        }
    }

    return lastError_;
}
