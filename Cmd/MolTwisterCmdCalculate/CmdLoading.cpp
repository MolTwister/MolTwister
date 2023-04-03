//
// Copyright (C) 2023 Richard Olsen.
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

#include "CmdLoading.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdLoading::getCmd()
{
    return "loading";
}

std::vector<std::string> CCmdLoading::getCmdLineKeywords()
{
    return { "loading" };
}

std::vector<std::string> CCmdLoading::getCmdHelpLines()
{
    return {
                "loading <DCD filename> <frame from> <frame to> <atom IDs to find loading for> <lower vector - loading region> <upper vector - loading region>"
           };
}

std::string CCmdLoading::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the loading curve from <frame from> to <frame to> of the atoms specified in\r\n";
    text+= "\t<atom IDs to find loading for> within the DCD file given by <DCD filename>. For each\r\n";
    text+= "\tframe, which is one point on the curve, the number of atoms from the list of atom IDs\r\n";
    text+= "\t(comma separated list without spaces, e.g., H,O,C7) that are within the loadin region\r\n";
    text+= "\tare counted. The loading region is specified by a lower and upper vector, which both\r\n";
    text+= "\tare entered as three numbers, <x> <y> <z>, separated by a space.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. <list of atoms for which loading is calculated>\r\n";
    text+= "\t2. <loading frame 1>\r\n";
    text+= "\t3. <loading frame 2>\r\n";
    text+= "\t      .\r\n";
    text+= "\t      .\r\n";
    text+= "\t      .\r\n";
    text+= "\tN+1. <loading frame N>\r\n";
    text+= "\twhere N is the number of included frames.";

    return text;
}

std::string CCmdLoading::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToInclude;
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

    if(frameTo <= frameFrom)
    {
        lastError_ = "Error: frame to must be larger than frame from!";
        return lastError_;
    }


    // Retrieve indices of loading targets (i.e. atoms that we wish to estimate loading of)
    std::string atomsToIncludeString;
    atomsToIncludeString = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(atomsToIncludeString);
    atomsToInclude = CASCIIUtility::getWords(atomsToIncludeString, ",");

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


    // Retrieve loading region
    C3DRect loadingRegion;
    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rLow_.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rLow_.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rLow_.z_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rHigh_.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rHigh_.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    loadingRegion.rHigh_.z_ = atof(text.data());


    // Prepare loading curve
    std::vector<int> loadingCurve;
    loadingCurve.resize(frameTo - frameFrom);


    // Go through all frames
    int index = 0;
    int numFrames = dcdFile.getNumRecords();
    pb.beginProgress("Calculating loading curves");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        int loading = 0;
        for(int i=0; i<atomIndicesToInclude.size(); i++)
        {
            C3DVector r = dcdFile.getCoordinate(atomIndicesToInclude[i]);
            if(loadingRegion.isWithin(r)) loading++;
        }

        loadingCurve[index++] = loading;

        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Print loading curve
    printf("\r\n");
    fprintf(stdOut_, "\t%-15s\r\n", atomsToIncludeString.data());
    for(int i=0; i<loadingCurve.size(); i++)
    {
        fprintf(stdOut_, "\t%-15i\r\n", loadingCurve[i]);
    }

    return lastError_;
}
