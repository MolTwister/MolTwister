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

#include "CmdBondLengthDyn.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdBondLengthDyn::getCmd()
{
    return "bondlengthdyn";
}

std::vector<std::string> CCmdBondLengthDyn::getCmdLineKeywords()
{
    return { "bondlengthdyn" };
}

std::vector<std::string> CCmdBondLengthDyn::getCmdHelpLines()
{
    return {
                "bondlengthdyn <DCD filename> <frame from> <frame to> <atom indices list (each pair defines a bond)>"
           };
}

std::string CCmdBondLengthDyn::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures distances between atom pairs defined in <atom indices list>, which is comma\r\n";
    text+= "\tseparated with no spaces alowed. E.g., <atom indices list>=0,3,5,6 will yield a\r\n";
    text+= "\tmeasurement of distances between atoms 0 and 3, as well as betwen atoms 5 and 6. More\r\n";
    text+= "\tdistances are added by extending the list. The distances are measured from the DCD\r\n";
    text+= "\tfile, <DCD filename>, for each frame between <frame from> and <frame to>. Thus, a\r\n";
    text+= "\tdistance plot, as function of time, is provided for each selected sistance to measure.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. <distance 1> <distance 2> ... <distance n>\r\n";
    text+= "\t2. <distance 1> <distance 2> ... <distance n>\r\n";
    text+= "\t         .\r\n";
    text+= "\t         .\r\n";
    text+= "\t         .\r\n";
    text+= "\tN. <distance 1> <distance 2> ... <distance n>\r\n";
    text+= "\twhere N is the number of selected frames from the DCD file and n is the number of pairs\r\n";
    text+= "\tof atoms in <atom indices list>, where the order is the same as in <atom indices list>.";

    return text;
}

std::string CCmdBondLengthDyn::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndices1;
    std::vector<int> atomIndices2;
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


    // Retrieve array of bonds we wish to measure the length of
    // Format: contigous list of atom indices, where each consequtive
    // pair of indices form a bond
    std::string stringAtomsToInclude;
    stringAtomsToInclude = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(stringAtomsToInclude);
    atomsToInclude = CASCIIUtility::getWords(stringAtomsToInclude, ",");

    for(int i=0; i<atomsToInclude.size(); i++)
    {
        if((i+1) >= atomsToInclude.size()) continue;

        atomIndices1.emplace_back((int)atoi(atomsToInclude[i++].data()));
        atomIndices2.emplace_back((int)atoi(atomsToInclude[i].data()));
    }


    // Print header
    char text1[64];
    if(stdOut_ == stdout) printf("\r\n");
    fprintf(stdOut_, "\t");
    for(int i=0; i<atomIndices1.size(); i++)
    {
        int iIndex1 = atomIndices1[i];
        int iIndex2 = atomIndices2[i];

        if(iIndex1 >= state_->atoms_.size()) continue;
        if(iIndex2 >= state_->atoms_.size()) continue;

        std::string ID1 = state_->atoms_[iIndex1]->getID();
        std::string ID2 = state_->atoms_[iIndex2]->getID();

        sprintf(text1, "%s-%s", ID1.data(), ID2.data());
        fprintf(stdOut_, "%-15s", text1);
    }
    fprintf(stdOut_, "\r\n");


    // Iterate through selected frames
    int numFrames = dcdFile.getNumRecords();
    if(stdOut_ != stdout) pb.beginProgress("Measure bond length dynamics");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        // Iterate through selected bonds to measure
        fprintf(stdOut_, "\t");
        for(int i=0; i<atomIndices1.size(); i++)
        {
            C3DVector r1 = dcdFile.getCoordinate(atomIndices1[i]);
            C3DVector r2 = dcdFile.getCoordinate(atomIndices2[i]);
            C3DVector r12 = r2 - r1;
            double R12 = r12.norm();

            fprintf(stdOut_, "%-15.3f", R12);
        }
        fprintf(stdOut_, "\r\n");

        if(stdOut_ != stdout) pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    if(stdOut_ != stdout) pb.endProgress();

    return lastError_;
}
