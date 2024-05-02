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

#include "CmdDihedralDistrCOM.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolecularTools.h"
#include "../Tools/MolecularSystemTools.h"

std::string CCmdDihedralDistrCOM::getCmd()
{
    return "dihedraldistrcom";
}

std::vector<std::string> CCmdDihedralDistrCOM::getCmdLineKeywords()
{
    return { "dihedraldistrcom", "indices", "atomtypes", "ignorepbc" };
}

std::vector<std::string> CCmdDihedralDistrCOM::getCmdHelpLines()
{
    return {
                "dihedraldistrcom <DCD filename> <frame from> <frame to> <num bins> <num COM bins> <COM dir> <COM range> indices <index 1> <index 2> <index 3> <index 4> [ignorepbc] indices <M> <COM mol index 1> ... <COM mol index M>",
                "dihedraldistrcom <DCD filename> <frame from> <frame to> <num bins> <num COM bins> <COM dir> <COM range> indices <index 1> <index 2> <index 3> <index 4> [ignorepbc] atomtypes <M> <COM mol atom ID 1> ... <COM mol atom ID M>",
                "dihedraldistrcom <DCD filename> <frame from> <frame to> <num bins> <num COM bins> <COM dir> <COM range> atomtypes <atom ID 1> <atom ID 2> <atom ID 3> <atom ID 4> <bond cutoff r> [ignorepbc] indices <M> <COM mol index 1> ... <COM mol index M>",
                "dihedraldistrcom <DCD filename> <frame from> <frame to> <num bins> <num COM bins> <COM dir> <COM range> atomtypes <atom ID 1> <atom ID 2> <atom ID 3> <atom ID 4> <bond cutoff r> [ignorepbc] atomtypes <M> <COM mol atom ID 1> ... <COM mol atom ID M>"
           };
}

std::string CCmdDihedralDistrCOM::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the dihedral distribution, based on the contents of the <DCD filename> DCD file,\r\n";
    text+= "\tfrom frame index <frame from> to frame index <frame to>. The number of bins to divide the full\r\n";
    text+= "\t360 degree angle into is defined by <num bins>. It is possible to either define four atom\r\n";
    text+= "\tindices (zero based) and hence study the conformational distributions of a single dihedral\r\n";
    text+= "\tover time, or to specify four atom IDs (e.g., H, O, C5) and study a specific type of dihedral\r\n";
    text+= "\tover time (i.e., several dihedrals of the same type, per frame). In the latter case, one also\r\n";
    text+= "\tneed to specify the bond cutoff radius, <bond cutoff r>, which defines how long a bond should\r\n";
    text+= "\tbe before it is no longer considered a bond. By default, bonds reach across periodic boundaries\r\n";
    text+= "\t(PBC), but it is possible to only consider bonds within the PBC by specifying 'ignorepbc'. \r\n";
    text+= "\r\n";
    text+= "\tIn the output from this calculation, the bins will be identified by the center of mass (COM)\r\n";
    text+= "\tpositions of the molecules containing the dihedral angles that are counted. Thus, the COM molecule\r\n";
    text+= "\tmust be identified, as well as the range of COM position that are to be binned and the number of\r\n";
    text+= "\tsuch bins. This is done through either specifying the 'indices' or 'atomtypes', together with the\r\n";
    text+= "\tnumber of indices or molecules, <M>, and the indices, <COM mol index i>, or atom IDs (such as H,\r\n";
    text+= "\tO, C6), <COM mol atom ID i>. The COM bins are defined by <num COM bins> <COM range> (of the form\r\n";
    text+= "\t<start pos> <end pos>), as well as the axis along which to perform binning <COM dir> in {x, y, z}.\r\n";
    text+= "\t\r\n";
    text+= "\tThe output from the calculation is given below.\r\n";
    text+= "\r\n";
    text+= "\tThe output starts with the following header:\r\n";
    text+= "\t1. Num angle Bins: <applied number of bins>\r\n";
    text+= "\t2. Num COM Bins: <applied number of bins>\r\n";
    text+= "\t3. Start COM range: <applied start pos>\r\n";
    text+= "\t4. End COM range: <applied end pos>\r\n";
    text+= "\t5. COM range direction (0=x,1=y,2=z): <applied direction>\r\n";
    text+= "\t6. From frame: <starting frame index>\r\n";
    text+= "\t7. To frame: <to frame index>\r\n";
    text+= "\t8. Bond cutoff: <applied bond cutoff>\r\n";
    text+= "\t9. Dihedral: <specified dihedral that was studied>\r\n";
    text+= "\t10. COMPos Angle Tot.Count Normalized\r\n";
    text+= "\r\n";
    text+= "\tFollowing the header is the data section:\r\n";
    text+= "\t11. <COM pos> <angle in degrees> <number of entries in this bin> <normalized num. entries, max is unity>\r\n";
    text+= "\t12. <COM pos> <angle in degrees> <number of entries in this bin> <normalized num. entries, max is unity>\r\n";
    text+= "\t                      .\r\n";
    text+= "\t                      .\r\n";
    text+= "\t                      .\r\n";
    text+= "\tNa*Nb+10. <COM pos> <angle in degrees> <number of entries in this bin> <normalized num. entries, max is unity>\r\n";
    text+= "\twhere Na is the number of COM bins and Nb is the number of angle bins (i.e., this constitutes.\r\n";
    text+= "\ta surface plot with angle along one axis and molecular COM along the other).\r\n";

    return text;
}

std::string CCmdDihedralDistrCOM::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToInclude[4];
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


    // Retrieve the number of bins to divide 360 degrees into
    text = CASCIIUtility::getArg(arguments, arg++);
    int numBins = atoi(text.data());


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
        lastError_= "Error: expected 'x', 'y' or 'z'!";
        return lastError_;
    }


    // Retrieve the COM range in previously chosen direction
    text = CASCIIUtility::getArg(arguments, arg++);
    double startRangeCOM = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double endRangeCOM = atof(text.data());
    if(fabs(startRangeCOM - endRangeCOM) < double(FLT_MIN))
    {
        lastError_ = "Error: cannot have a start range equal to end range!";
        return lastError_;
    }


    // Retrieve atom types for which to calculate dihedral angle distributions
    C3DRect pbc;
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
    pbc = state_->view3D_->calcPBC(state_->currentFrame_);
    bool distAcrossPBC = true;

    std::string arrayOfIDs[4];
    double bondCutoffSqr = -1.0;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "indices")
    {
        arrayOfIDs[0] = CASCIIUtility::getArg(arguments, arg++);
        atomIndicesToInclude[0].emplace_back(atoi(arrayOfIDs[0].data()));

        arrayOfIDs[1] = CASCIIUtility::getArg(arguments, arg++);
        atomIndicesToInclude[1].emplace_back(atoi(arrayOfIDs[1].data()));

        arrayOfIDs[2] = CASCIIUtility::getArg(arguments, arg++);
        atomIndicesToInclude[2].emplace_back(atoi(arrayOfIDs[2].data()));

        arrayOfIDs[3] = CASCIIUtility::getArg(arguments, arg++);
        atomIndicesToInclude[3].emplace_back(atoi(arrayOfIDs[3].data()));

        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "ignorepbc") distAcrossPBC = false;
        else arg--;
    }
    else if(text == "atomtypes")
    {
        std::string ID;
        std::map<CAtom*,int> mapAtomPtrToIndex;
        std::vector<CAtom*> atomsType1, atomsType2, atomsType3, atomsType4;

        state_->getMapAtomPtrToIndex(mapAtomPtrToIndex);

        arrayOfIDs[0] = CASCIIUtility::getArg(arguments, arg++);
        arrayOfIDs[1] = CASCIIUtility::getArg(arguments, arg++);
        arrayOfIDs[2] = CASCIIUtility::getArg(arguments, arg++);
        arrayOfIDs[3] = CASCIIUtility::getArg(arguments, arg++);

        state_->getAtomsWithID(arrayOfIDs[0], atomsType1);
        state_->getAtomsWithID(arrayOfIDs[1], atomsType2);
        state_->getAtomsWithID(arrayOfIDs[2], atomsType3);
        state_->getAtomsWithID(arrayOfIDs[3], atomsType4);

        text = CASCIIUtility::getArg(arguments, arg++);
        bondCutoffSqr = atof(text.data())*atof(text.data());

        text = CASCIIUtility::getArg(arguments, arg++);
        if(text == "ignorepbc") distAcrossPBC = false;
        else arg--;

        pb.beginProgress("Searching for dihedrals");

        CMolecularTools molTools(state_, stdOut_);
        for(int j=0; j<atomsType1.size(); j++)
        {
            for(int k=0; k<atomsType2.size(); k++)
            {
                if(atomsType1[j] == atomsType2[k]) continue;
                if(molTools.calcDistance2(atomsType1[j], atomsType2[k], state_->currentFrame_, pbc, distAcrossPBC) < bondCutoffSqr)
                {
                    for(int l=0; l<atomsType3.size(); l++)
                    {
                        if(atomsType3[l] == atomsType2[k]) continue;
                        if(atomsType3[l] == atomsType1[j]) continue;
                        if(molTools.calcDistance2(atomsType2[k], atomsType3[l], state_->currentFrame_, pbc, distAcrossPBC) < bondCutoffSqr)
                        {
                            for(int m=0; m<atomsType4.size(); m++)
                            {
                                if(atomsType4[m] == atomsType3[l]) continue;
                                if(atomsType4[m] == atomsType2[k]) continue;
                                if(atomsType4[m] == atomsType1[j]) continue;
                                if(molTools.calcDistance2(atomsType3[l], atomsType4[m], state_->currentFrame_, pbc, distAcrossPBC) < bondCutoffSqr)
                                {
                                    int index1 = mapAtomPtrToIndex[atomsType1[j]];
                                    int index2 = mapAtomPtrToIndex[atomsType2[k]];
                                    int index3 = mapAtomPtrToIndex[atomsType3[l]];
                                    int index4 = mapAtomPtrToIndex[atomsType4[m]];

                                    atomIndicesToInclude[0].emplace_back(index1);
                                    atomIndicesToInclude[1].emplace_back(index2);
                                    atomIndicesToInclude[2].emplace_back(index3);
                                    atomIndicesToInclude[3].emplace_back(index4);
                                }
                            }
                        }
                    }
                }
            }

            pb.updateProgress(j, (int)atomsType1.size());
        }

        pb.endProgress();

        CMolecularSystemTools(state_, stdOut_).removeDuplicateDihIndices(atomIndicesToInclude[0], atomIndicesToInclude[1], atomIndicesToInclude[2], atomIndicesToInclude[3], nullptr);
        printf("\r\n\tFound %i dihedrals of type %s-%s-%s-%s!\r\n", (int)atomIndicesToInclude[0].size(), arrayOfIDs[0].data(), arrayOfIDs[1].data(), arrayOfIDs[2].data(), arrayOfIDs[3].data());
        if(!distAcrossPBC) printf("\tNote that bonds across PBC were not included in search!\r\n");
    }
    else
    {
        lastError_ = "Error: expected 'indices' or 'atomtypes'!";
        return lastError_;
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


    // Go through DCD and calculate dihedral distributions of the selected dihedrals
    C3DRect* pbcPtr = nullptr;
    if(distAcrossPBC) pbcPtr = &pbc;
    int numProcessedFrames = 0;
    int numFrames = dcdFile.getNumRecords();
    std::vector<std::vector<int>> comAngleBins;
    std::vector<std::vector<double>> comAngleBinsNormalized;
    comAngleBins.resize(numCOMBins);
    comAngleBinsNormalized.resize(numCOMBins);
    for(int i=0; i<comAngleBins.size(); i++)
    {
        comAngleBins[i].resize(numBins, 0);
        comAngleBinsNormalized[i].resize(numBins, 0.0);
    }

    pb.beginProgress("Calculating dihedral-COM distributions");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        int index1;
        int index2;
        int index3;
        int index4;

        dcdFile.gotoRecord(t);
        CDCDTools dcdTools(state_, stdOut_);
        C3DVector vecCOM = dcdTools.getCenterOfMass(atomIndicesToIncludeCOM, &dcdFile);
        double com = 0.0;
        if(comDir == 0) com = vecCOM.x_;
        if(comDir == 1) com = vecCOM.y_;
        if(comDir == 2) com = vecCOM.z_;

        for(int i=0; i<atomIndicesToInclude[0].size();  i++)
        {
            index1 = atomIndicesToInclude[0][i];
            index2 = atomIndicesToInclude[1][i];
            index3 = atomIndicesToInclude[2][i];
            index4 = atomIndicesToInclude[3][i];

            double angle = dcdTools.measureDihedral(dcdFile, index1, index2, index3, index4, pbcPtr);
            int binIndex = int((angle / (2.0*M_PI)) * double(numBins));
            if(binIndex == numBins) binIndex = 0;

            int binCOMIndex = int((com - startRangeCOM) / (endRangeCOM - startRangeCOM) * double(numCOMBins));

            if((binIndex < numBins) && (binCOMIndex >= 0) && (binCOMIndex < numCOMBins)) comAngleBins[binCOMIndex][binIndex]++;
        }

        numProcessedFrames++;
        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Normalize the dihedral distribution in each bin, where
    // maximum is set to unity
    for(int i=0; i<comAngleBins.size(); i++)
    {
        int max = 0;

        // Find maximum
        for(int j=0; j<comAngleBins[i].size(); j++)
        {
            if(comAngleBins[i][j] > max) max = comAngleBins[i][j];
        }

        // Create normalized value
        for(int j=0; j<comAngleBins[i].size(); j++)
        {
            if(max > 0) comAngleBinsNormalized[i][j] = double(comAngleBins[i][j]) / double(max);
            else        comAngleBinsNormalized[i][j] = 0.0;
        }
    }


    // Present the results
    printf("\r\n");
    fprintf(stdOut_, "\tNum angle bins: %i\r\n", numBins);
    fprintf(stdOut_, "\tNum COM bins: %i\r\n", numCOMBins);
    fprintf(stdOut_, "\tStart COM range: %g\r\n", startRangeCOM);
    fprintf(stdOut_, "\tEnd COM range: %g\r\n", endRangeCOM);
    fprintf(stdOut_, "\tCOM range direction (0=x,1=y,2=z): %i\r\n", comDir);
    fprintf(stdOut_, "\tFrom frame: %i\r\n", frameFrom);
    fprintf(stdOut_, "\tTo frame: %i\r\n", frameTo);
    fprintf(stdOut_, "\tBond cutoff: %.6f\r\n", sqrt(bondCutoffSqr));
    fprintf(stdOut_, "\tDihedral: %s-%s-%s-%s\r\n\r\n", arrayOfIDs[0].data(), arrayOfIDs[1].data(), arrayOfIDs[2].data(), arrayOfIDs[3].data());
    fprintf(stdOut_, "\t%-15s%-15s%-15s%-15s\r\n", "COMPos", "Angle", "Tot.Count", "Normalized");
    for(int i=0; i<numCOMBins; i++)
    {
        for(int j=0; j<numBins; j++)
        {
            fprintf(stdOut_, "\t% -15.4g% -15.4g%-15i%-15.3f\r\n", 0.5 * (endRangeCOM - startRangeCOM) / double(numCOMBins) * (2.0*double(i) + 1.0) + startRangeCOM, 0.5 * (2.0*M_PI / double(numBins)) * (2.0*double(j) + 1.0) * 180.0 / M_PI, comAngleBins[i][j], comAngleBinsNormalized[i][j]);
        }
    }

    return lastError_;
}
