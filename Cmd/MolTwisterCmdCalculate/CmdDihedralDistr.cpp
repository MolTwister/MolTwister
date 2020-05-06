#include "CmdDihedralDistr.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolecularTools.h"
#include "../Tools/MolecularSystemTools.h"

std::string CCmdDihedralDistr::getCmd()
{
    return "dihedraldistr";
}

std::vector<std::string> CCmdDihedralDistr::getCmdLineKeywords()
{
    return { "dihedraldistr", "indices", "atomtypes", "ignorepbc" };
}

std::vector<std::string> CCmdDihedralDistr::getCmdHelpLines()
{
    return {
                "dihedraldistr <DCD filename> <frame from> <frame to> <num bins> indices <index 1> <index 2> <index 3> <index 4>",
                "dihedraldistr <DCD filename> <frame from> <frame to> <num bins> atomtypes <atom ID 1> <atom ID 2> <atom ID 3> <atom ID 4> <bond cutoff r> [ignorepbc] "
           };
}

std::string CCmdDihedralDistr::getCmdFreetextHelp()
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
    text+= "\t(PBC), but it is possible to only consider bonds within the PBC by specifying 'ignorepbc'. The\r\n";
    text+= "\toutput from the calculation is given below.\r\n";
    text+= "\r\n";
    text+= "\tThe output starts with the following header:\r\n";
    text+= "\t1. Num Bins: <applied number of bins>\r\n";
    text+= "\t2. From frame: <starting frame index>\r\n";
    text+= "\t3. To frame: <to frame index>\r\n";
    text+= "\t4. Bond cutoff: <applied bond cutoff>\r\n";
    text+= "\t5. Dihedral: <specified dihedral that was studied>\r\n";
    text+= "\t6. Index Angle Tot.Count Count/Frm\r\n";
    text+= "\r\n";
    text+= "\tFollowing the header is the data section:\r\n";
    text+= "\t7. <bin index> <angle in degrees> <number of entries in this bin> <number of entries per frame in this bin>\r\n";
    text+= "\t8. <bin index> <angle in degrees> <number of entries in this bin> <number of entries per frame in this bin>\r\n";
    text+= "\t                      .\r\n";
    text+= "\t                      .\r\n";
    text+= "\t                      .\r\n";
    text+= "\tN+6. <bin index> <angle in degrees> <number of entries in this bin> <number of entries per frame in this bin>\r\n";
    text+= "\twhere N is the number of bins.";

    return text;
}

std::string CCmdDihedralDistr::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<int> atomIndicesToInclude[4];
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


    // Retrieve atom types for which to calculate dihedral angle distributions
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
    }
    else if(text == "atomtypes")
    {
        bool distAcrossPBC = true;
        C3DRect pbc;
        std::string ID;
        std::map<CAtom*,int> mapAtomPtrToIndex;
        std::vector<CAtom*> atomsType1, atomsType2, atomsType3, atomsType4;

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


    // Go through DCD and calculate dihedral distributions of the selected dihedrals
    int numProcessedFrames = 0;
    int numFrames = dcdFile.getNumRecords();
    std::vector<int> angleBins;
    angleBins.resize(numBins, 0);
    pb.beginProgress("Calculating dihedral distributions");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            printf("Error: frame %i does not exist (num. frames = %i)", t, numFrames);
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        int index1;
        int index2;
        int index3;
        int index4;

        dcdFile.gotoRecord(t);
        for(int i=0; i<atomIndicesToInclude[0].size();  i++)
        {
            index1 = atomIndicesToInclude[0][i];
            index2 = atomIndicesToInclude[1][i];
            index3 = atomIndicesToInclude[2][i];
            index4 = atomIndicesToInclude[3][i];

            double angle = CDCDTools(state_, stdOut_).measureDihedral(dcdFile, index1, index2, index3, index4);
            int bin = int((angle / (2.0*M_PI)) * double(numBins));
            if(bin == numBins) bin = 0;
            if(bin < numBins) angleBins[bin]++;
        }

        numProcessedFrames++;
        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Calculate the number of angles found per frame
    std::vector<double> angleBinsPerFrame;
    angleBinsPerFrame.resize(numBins, 0.0);
    for(int i=0; i<numBins; i++)
    {
        if(numProcessedFrames != 0)
            angleBinsPerFrame[i] = angleBins[i] / double(numProcessedFrames);
    }


    // Present the results
    printf("\r\n");
    fprintf(stdOut_, "\tNum Bins: %i\r\n", numBins);
    fprintf(stdOut_, "\tFrom frame: %i\r\n", frameFrom);
    fprintf(stdOut_, "\tTo frame: %i\r\n", frameTo);
    fprintf(stdOut_, "\tBond cutoff: %.6f\r\n", sqrt(bondCutoffSqr));
    fprintf(stdOut_, "\tDihedral: %s-%s-%s-%s\r\n\r\n", arrayOfIDs[0].data(), arrayOfIDs[1].data(), arrayOfIDs[2].data(), arrayOfIDs[3].data());
    fprintf(stdOut_, "\t%-15s%-15s%-15s%15s\r\n", "Index", "Angle", "Tot.Count", "Count/Frm");
    for(int i=0; i<numBins; i++)
    {
        fprintf(stdOut_, "\t%-15i%-15.4f%-15i%-15.8f\r\n", i, 0.5 * (2.0*M_PI / double(numBins)) * (2.0*double(i) + 1.0) * 180.0 / M_PI, angleBins[i], angleBinsPerFrame[i]);
    }

    return lastError_;
}
