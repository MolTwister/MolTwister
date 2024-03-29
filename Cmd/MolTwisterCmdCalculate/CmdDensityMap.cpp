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

#include "CmdDensityMap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdDensityMap::getCmd()
{
    return "densitymap";
}

std::vector<std::string> CCmdDensityMap::getCmdLineKeywords()
{
    return { "densitymap" };
}

std::vector<std::string> CCmdDensityMap::getCmdHelpLines()
{
    return {
                "densitymap <DCD filename> <stride> <span> <last frame to load> <num bins> <atoms (comma sep, no space)> <direction>"
           };
}

std::string CCmdDensityMap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates density maps of a specified atomic selection. Each density map that is generated will be\r\n";
    text+= "\ta count taken over <span> frames of the specified DCD file, centered at an index, t. The index, t,\r\n";
    text+= "\truns over all frames of the DCD file, up to <last frame to load>, where each index is separated by <stride>\r\n";
    text+= "\tframes. I.e., the number of generated density maps are approximately <last frame to load> / <stride>.\r\n";
    text+= "\tEach map is generated by specifying a plane, given by <direction> in {xy, yz, zx}, for which all atoms\r\n";
    text+= "\tare projected into. The atoms specified by the comma separated list of atom IDs (i.e., <atoms>, s.a.,\r\n";
    text+= "\tO, H, C5, etc.) are counted within every bin of the defined plane, where each side of the plane are\r\n";
    text+= "\tdivided into <num bins> equisized bins. The output is formatted as described below.\r\n";
    text+= "\r\n";
    text+= "\tFirst a header, consisting of three lines are output:\r\n";
    text+= "\t1. x(t,n,val) y(t,n,val) xbin(t,n,val) ybin(t,n,val) bincnt(t,n,val) ...\r\n";
    text+= "\t2. <t>        <t>        <t>           <t>           <t>             ...\r\n";
    text+= "\t3. <Nx>       <Ny>       <num bins>^2  <num bins>^2  <num bins>^2    ...\r\n";
    text+= "\tThe first line describes each column that will follow in the data section, located immediately below the\r\n";
    text+= "\theader. The ellipsis (...) denotes that the columns are repeated for each output frame, <t>, which is\r\n";
    text+= "\tcounted over <span>/2 on each side. Note that if the span is outside the available frames of the DCD\r\n";
    text+= "\tfile, then these frames are simply ignored. The number of points (i.e., the length of the data columns\r\n";
    text+= "\tbelow) are given by <Nx> and <Ny>., while the total number of bins in the specified plane is given by the\r\n";
    text+= "\tsquare of <num bins>. The data section is formatted as:\r\n";
    text+= "\t4. <x> <y> <x center pos. of bin containing (x,y)> <y center pos. of bin containing (x,y)> <bin count of bin containing (x,y) ...\r\n";
    text+= "\t5. <x> <y> <x center pos. of bin containing (x,y)> <y center pos. of bin containing (x,y)> <bin count of bin containing (x,y) ...\r\n";
    text+= "                                    .\r\n";
    text+= "                                    .\r\n";
    text+= "                                    .\r\n";
    text+= "\tN+3. <x> <y> <x center pos. of bin containing (x,y)> <y center pos. of bin containing (x,y)> <bin count of bin containing (x,y) ...\r\n";
    text+= "\twhere N is the number of selected atoms.";

    return text;
}

std::string CCmdDensityMap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToInclude;
    CDCDFile dcdFile;
    std::string text;
    int every, span, nBins, lastFrame;


    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);

    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
        return lastError_;
    }

    int numFrames = dcdFile.getNumRecords();
    text = CASCIIUtility::getArg(arguments, arg++);
    every = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    span = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    lastFrame = atoi(text.data());
    if(lastFrame == -1) lastFrame = numFrames;

    text = CASCIIUtility::getArg(arguments, arg++);
    nBins = atoi(text.data());
    int nBins2 = nBins*nBins;

    if(every <= 0)
    {
        lastError_ = "Error: stride must be larger than zero!";
        return lastError_;
    }
    if(nBins <= 0)
    {
        lastError_ = "Error: number of bins must be larger than zero!";
        return lastError_;
    }


    // Retrieve indices of density targets (i.e. atoms that we wish to estimate density of)
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


    // Retrieve direction
    int direction = 0;
    std::string directionString;
    directionString = CASCIIUtility::getArg(arguments, arg++);
    if(directionString == "xy")          direction = 0;
    else if(directionString == "yz")     direction = 1;
    else if(directionString == "zx")     direction = 2;
    else
    {
        lastError_ = std::string("Error: expected 'x', 'y', or 'z', not ") + directionString + std::string("!");
        return lastError_;
    }


    // Prepare arrays to collect points, as well as bin densities
    int numDensityMaps = (int)ceil(double(lastFrame) / double(every));
    std::vector<std::vector<double>> pointArraysX;
    std::vector<std::vector<double>> pointArraysY;
    std::vector<std::vector<int>> aBinCount;
    pointArraysX.resize(numDensityMaps);
    pointArraysY.resize(numDensityMaps);
    aBinCount.resize(numDensityMaps);
    for(int i=0; i<aBinCount.size(); i++)
    {
        aBinCount[i].resize(nBins2, 0);
    }


    // Calculate parameters used in binning process
    C3DRect pbc = state_->view3D_->calcPBC();
    double rangeX, rangeY, XMin, YMin;
    if(direction == 0)         { XMin = pbc.rLow_.x_; YMin = pbc.rLow_.y_; rangeX = pbc.getWidthX(); rangeY = pbc.getWidthY(); }
    else if(direction == 1)    { XMin = pbc.rLow_.y_; YMin = pbc.rLow_.z_; rangeX = pbc.getWidthY(); rangeY = pbc.getWidthZ(); }
    else if(direction == 2)    { XMin = pbc.rLow_.z_; YMin = pbc.rLow_.x_; rangeX = pbc.getWidthZ(); rangeY = pbc.getWidthX(); }
    else                       { XMin = pbc.rLow_.z_; YMin = pbc.rLow_.x_; rangeX = pbc.getWidthZ(); rangeY = pbc.getWidthX(); }

    if((rangeX == 0.0) || (rangeY == 0.0))
    {
        lastError_ = std::string("Error: found zero PBC range (") + std::to_string(rangeX) + std::string(", ") + std::to_string(rangeY) + std::string(")!");
        return lastError_;
    }

    double binCoeffX = double(nBins) / rangeX;
    double binCoeffY = double(nBins) / rangeY;


    // Go through all frames
    int span_2 =span / 2;
    int index = 0;
    pb.beginProgress("Calculate densitymaps");
    for(int i=span_2; i<lastFrame; i+=every)
    {
        for(int t=(i-span_2); t<(i+span_2); t++)
        {
            if((t < 0) || (t >= numFrames)) continue;
            dcdFile.gotoRecord(t);

            for(int j=0; j<atomIndicesToInclude.size(); j++)
            {
                C3DVector r = dcdFile.getCoordinate(atomIndicesToInclude[j]);

                double pos[2];
                if(direction == 0)         { pos[0] = r.x_; pos[1] = r.y_; } // xy
                else if(direction == 1)    { pos[0] = r.y_; pos[1] = r.z_; } // yz
                else if(direction == 2)    { pos[0] = r.z_; pos[1] = r.x_; } // zx
                else                       { pos[0] = r.z_; pos[1] = r.x_; } // zx

                pointArraysX[index].emplace_back(pos[0]);
                pointArraysY[index].emplace_back(pos[1]);

                int XBin = (int)(binCoeffX * (pos[0] - XMin));
                int YBin = (int)(binCoeffY * (pos[1] - YMin));

                if((XBin >= 0) && (XBin <= nBins) && (YBin >= 0) && (YBin <= nBins))
                {
                    int bin = YBin*nBins + XBin;
                    aBinCount[index][bin]++;
                }
            }
        }

        index++;
        pb.updateProgress(i, numFrames);
    }
    pb.endProgress();


    // Find lenght of longest column to be printed
    int lengthOfLongestColumn = -1;
    for(int i=0; i<numDensityMaps; i++)
    {
        if(int(pointArraysX[i].size()) > lengthOfLongestColumn) lengthOfLongestColumn = (int)pointArraysX[i].size();
        if(int(pointArraysY[i].size()) > lengthOfLongestColumn) lengthOfLongestColumn = (int)pointArraysY[i].size();
        if(int(aBinCount[i].size()) > lengthOfLongestColumn) lengthOfLongestColumn = (int)aBinCount[i].size();
    }

    if(lengthOfLongestColumn == -1)
    {
        lastError_ = "Error: no density profile generated!";
        return lastError_;
    }


    // Print header of results (type, frame and length)
    printf("\r\n");
    fprintf(stdOut_, "\t");
    for(int i=0; i<numDensityMaps; i++) { fprintf(stdOut_, "%-15s%-15s%-15s%-15s%-17s", "x(t,n,val)", "y(t,n,val)", "xbin(t,n,val)", "ybin(t,n,val)", "bincnt(t,n,val)"); }

    fprintf(stdOut_, "\r\n\t");
    for(int i=0; i<numDensityMaps; i++)
    {
        int t = span_2 + i*every;
        fprintf(stdOut_, "%-15i%-15i%-15i%-15i%-17i", t, t, t, t, t);
    }

    fprintf(stdOut_, "\r\n\t");
    for(int i=0; i<numDensityMaps; i++)
    {
        int iNx = (int)pointArraysX[i].size();
        int iNy = (int)pointArraysY[i].size();
        fprintf(stdOut_, "%-15i%-15i%-15i%-15i%-17i", iNx, iNy, nBins2, nBins2, nBins2);
    }


    // Print results
    double deltaX = rangeX / double(nBins);
    double deltaY = rangeY / double(nBins);
    fprintf(stdOut_, "\r\n");
    for(int i=0; i<lengthOfLongestColumn; i++)
    {
        fprintf(stdOut_, "\t");
        for(int j=0; j<numDensityMaps; j++)
        {
            int Nx = (int)pointArraysX[j].size();
            int Ny = (int)pointArraysY[j].size();

            if((i >= Nx) || (i >= Ny)) fprintf(stdOut_, "% -15g% -15g", 0.0, 0.0);
            else
            {
                fprintf(stdOut_, "% -15g% -15g", pointArraysX[j][i], pointArraysY[j][i]);
            }

            if(i >= nBins2) fprintf(stdOut_, "% -15g% -15g% -17i", 0.0, 0.0, 0);
            else
            {
                // Since i = iYBin*iNBins + iXBin, we have that
                // iY = (i div iNBins), and that iX = (i mod iNBins)
                int Y = i / nBins;
                int X = i % nBins;
                fprintf(stdOut_, "% -15g% -15g% -17i", XMin + (double(X)+0.5)*deltaX, YMin + (double(Y)+0.5)*deltaY, aBinCount[j][i]);
            }
        }
        fprintf(stdOut_, "\r\n");
    }
    return lastError_;
}
