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

#include "CmdHBondCount.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdHBondCount::getCmd()
{
    return "hbondcount";
}

std::vector<std::string> CCmdHBondCount::getCmdLineKeywords()
{
    return { "hbondcount", "stride", "pbcfromvisual", "searchtype", "intermolecular", "intramolecular", "allhbonds", "nopbc" };
}

std::vector<std::string> CCmdHBondCount::getCmdHelpLines()
{
    return {
                "hbondcount <DCD filename> <frame from> <frame to> stride <stride> <M> <h-bond crit 1> ... <h-bond crit M> [pbcfromvisual] [searchtype <type>] [nopbc]",
                "hbondcount <DCD filename> <frame from> <frame to> <M> <h-bond crit 1> ... <h-bond crit M> [pbcfromvisual] [searchtype <type>] [nopbc]"
           };
}

std::string CCmdHBondCount::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates hydrogen bond statistics from frames <frame from> to <frame to> within the DCD file specified\r\n";
    text+= "\tby <DCD filenmae>. <M> hydrogen bond criteria are specified in the <h-bond crit i> parameters. Each such\r\n";
    text+= "\tparameter is given as a comma separated list with no spaces in-between. There are 7 entries in the comma\r\n";
    text+= "\tseparated lists\r\n";
    text+= "\t* the donor atom to consider (e.g., O, O4, C7)\r\n";
    text+= "\t* the hydrogen atom to consider\r\n";
    text+= "\t* the acceptor atom to consider\r\n";
    text+= "\t* the distance criteria between donor and hydrogen\r\n";
    text+= "\t* the distance criteria between hydrogen and acceptor\r\n";
    text+= "\t* the distance criteria between donor and acceptor\r\n";
    text+= "\t* the angle criteria, donor-hydrogen-acceptor\r\n";
    text+= "\tBoth distance and angular criteria can be specified as '-', which means they are ignored. It is possible to\r\n";
    text+= "\tload every <stride> frame from the DCD file, which can be useful for large DCD files. By default, the periodic\r\n";
    text+= "\tboundary conditions (PBC) are taken from each frame in the DCD file. However, by specifying 'pbcfromvisual'\r\n";
    text+= "\tthe PBC is collected from the one active in the 3D view. Using 'nopbc' ignores PBCs. By default, the bonds\r\n";
    text+= "\tare searched by considereing only intermolecular hydrogen bonds. By using the 'searchtype' parameter, it is\r\n";
    text+= "\tpossible to specify\r\n";
    text+= "\t* <type>=intermolecular: only consider intermolecular bonds\r\n";
    text+= "\t* <type>=intramolecular: only consider intramolecular bonds\r\n";
    text+= "\t* <type>=allhbonds: consider all types of bonds\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. [1] D:<donor> H:<hydrogen> A:<acceptor> d-h:<dist.crit. donor-hydrogen> h-a:<dist.crit. hydrogen-acceptor> d-a:<dist.crit. donor-accepror> d-h...a:<angle crit>\r\n";
    text+= "\t2. [2] D:<donor> H:<hydrogen> A:<acceptor> d-h:<dist.crit. donor-hydrogen> h-a:<dist.crit. hydrogen-acceptor> d-a:<dist.crit. donor-accepror> d-h...a:<angle crit>\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\tN. [N] D:<donor> H:<hydrogen> A:<acceptor> d-h:<dist.crit. donor-hydrogen> h-a:<dist.crit. hydrogen-acceptor> d-a:<dist.crit. donor-accepror> d-h...a:<angle crit>\r\n";
    text+= "\tN+1. Index Tot.count Num.Donors Num.Accept. Num.Conn.Don. Num.Conn.Acc. Frame FirstDon. FirstHydr. FirstAcc\r\n";
    text+= "\tN+2. <index> <tot.count> <num. donors> <num. acceptors> <num h-bonds conn to donors> <num h-bonds conn to acceptors> <frame index> <first donor in frame> <first hydrogen in frame> <first acceptor in frame>\r\n";
    text+= "\tN+3. <index> <tot.count> <num. donors> <num. acceptors> <num h-bonds conn to donors> <num h-bonds conn to acceptors> <frame index> <first donor in frame> <first hydrogen in frame> <first acceptor in frame>\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\t            .\r\n";
    text+= "\tN+K+1. <index> <tot.count> <num. donors> <num. acceptors> <num h-bonds conn to donors> <num h-bonds conn to acceptors> <frame index> <first donor in frame> <first hydrogen in frame> <first acceptor in frame>\r\n";
    text+= "\twhere N is the number of hydrogen bond criteria and K is the number of analyzed frames.";

    return text;
}

std::string CCmdHBondCount::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    CDCDFile dcdFile;
    std::string text;
    int frameFrom, frameTo, stride;


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

    stride = 1;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "stride")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        stride = atoi(text.data());

    } else arg--;

    if(stride < 1)
    {
        lastError_ = "Error: stride cannot be less than 1!";
        return lastError_;
    }

    if((frameTo - frameFrom) <= 0)
    {
        lastError_ = "Error: no frames were selected!";
        return lastError_;
    }


    // Set up H-bond criteria (num criteria, comma sep list of criteria 1, ....)
    // Comma sep list is without space and consists of D, H, A, r_dh_c, r_ha_c, r_da_c, alpha_dha_c,
    // where r and alpha may also be '-', which means "ignore" and leads to DBL_MAX values (and 0.0 for angle)
    std::vector<CDCDTools::CHBondCriteria> HBondCriteria;

    text = CASCIIUtility::getArg(arguments, arg++);
    int numCriteria = atoi(text.data());
    for(int i=0; i<numCriteria; i++)
    {
        std::vector<std::string> criteria;
        std::string criteriaString;
        criteriaString = CASCIIUtility::getArg(arguments, arg++);
        CASCIIUtility::removeWhiteSpace(criteriaString);
        criteria = CASCIIUtility::getWords(criteriaString, ",");

        if(criteria.size() != 7)
        {
            lastError_ = "Error: each set of H-bond criteria must consist of 7 numbers (D, H, A, r_dh_c, r_ha_c, r_da_c, alpha_dha_c)!";
            return lastError_;
        }

        CDCDTools::CHBondCriteria crit;
        crit.setDonor(criteria[0]);
        crit.setHydrogen(criteria[1]);
        crit.setAcceptor(criteria[2]);

        if(criteria[3] == "-")   crit.setR_dh(DBL_MAX);
        else                     crit.setR_dh(atof(criteria[3].data()));

        if(criteria[4] == "-")   crit.setR_ha(DBL_MAX);
        else                     crit.setR_ha(atof(criteria[4].data()));

        if(criteria[5] == "-")   crit.setR_da(DBL_MAX);
        else                     crit.setR_da(atof(criteria[5].data()));

        if(criteria[6] == "-")   crit.setAlpha_dha(0.0);
        else                     crit.setAlpha_dha(atof(criteria[6].data()));

        HBondCriteria.emplace_back(crit);
    }


    // Determine type of PBC
    bool usePBCFromVisual = false;
    C3DRect pbc;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "pbcfromvisual")
    {
        pbc = state_->view3D_->calcPBC();
        usePBCFromVisual = true;

    } else arg--;


    // Determine span of Hbonds (inter-molecular, intra-molecular or both)
    CDCDTools::EHBondSpan HBondSpan = CDCDTools::spanOnlyInterMolecularHBonds;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "searchtype")
    {
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "intermolecular") HBondSpan = CDCDTools::spanOnlyInterMolecularHBonds;
        else if(text == "intramolecular") HBondSpan = CDCDTools::spanOnlyIntraMolecularHBonds;
        else if(text == "allhbonds") HBondSpan = CDCDTools::spanAllHBonds;
        else
        {
            lastError_ = "Error: expected 'intermolecular', 'intramolecular', or 'allhbonds'!";
            return lastError_;
        }

    } else arg--;


    // Determine if we want to ignore PBCs
    bool noPBC = false;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "nopbc")
    {
        noPBC = true;

    } else arg--;


    // Loop through selected frames
    std::vector<int> HBondTotCounts;
    std::vector<int> numDonors;
    std::vector<int> numAcceptors;
    std::vector<int> numHBondsConnToDonorsArray;
    std::vector<int> numHBondsConnToAcceptorsArray;
    std::vector<CDCDTools::CHBond> firstHBondInFrame;
    int numFrames = dcdFile.getNumRecords();
    pb.beginProgress("Searching for H-bonds");
    for(int t=frameFrom; t<frameTo; t+=stride)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        if(!usePBCFromVisual)
        {
            pbc = dcdFile.getCurrentPBC();
        }


        // Find total number of H-bonds within the system
        std::shared_ptr<std::vector<CDCDTools::CHBond>> HBonds = CDCDTools(state_, stdOut_).retrieveHBondsCurrDCDFrame(&dcdFile, &pbc, HBondCriteria, HBondSpan, noPBC);
        HBondTotCounts.emplace_back((int)HBonds->size());


        // Store first HBond in list of H-bonds
        CDCDTools::CHBond firstHBond;
        if(HBonds->size() > 0) firstHBond = (*HBonds)[0];
        firstHBondInFrame.emplace_back(firstHBond);


        // Find total number of H-bond ends connected to donors (including donor hydrogens forming hydrogen bonds)
        // Start by finding all possible donor indices (but do not count them twice). Then, go through all donors
        // and count the number of H-bonds connected to each (and sum up all)
        std::vector<int> donors;
        for(int i=0; i<HBonds->size(); i++)
        {
            bool exists = false;
            for(int j=0; j<donors.size(); j++) { if(donors[j] == (*HBonds)[i].d_) exists = true; }
            if(!exists) donors.emplace_back((*HBonds)[i].d_);
        }
        numDonors.emplace_back((int)donors.size());

        int numHBondsConnToDonors = 0;
        for(int i=0; i<donors.size(); i++)
        {
            for(int j=0; j<HBonds->size(); j++)
            {
                if((*HBonds)[j].d_ == donors[i]) numHBondsConnToDonors++;
                if((*HBonds)[j].a_ == donors[i]) numHBondsConnToDonors++;
            }
        }
        numHBondsConnToDonorsArray.emplace_back(numHBondsConnToDonors);


        // Find total number of H-bond ends connected to acceptors (including donor hydrogens forming  hydrogen bonds)
        // Start by finding all possible acceptor indices (but do not count them twice). Then, go through all acceptors
        // and count the number of H-bonds connected to each (and sum up all)
        std::vector<int> acceptors;
        for(int i=0; i<HBonds->size(); i++)
        {
            bool exists = false;
            for(int j=0; j<acceptors.size(); j++) { if(acceptors[j] == (*HBonds)[i].a_) exists = true; }
            if(!exists) acceptors.emplace_back((*HBonds)[i].a_);
        }
        numAcceptors.emplace_back((int)acceptors.size());

        int numHBondsConnToAcceptors = 0;
        for(int i=0; i<acceptors.size(); i++)
        {
            for(int j=0; j<HBonds->size(); j++)
            {
                if((*HBonds)[j].d_ == acceptors[i]) numHBondsConnToAcceptors++;
                if((*HBonds)[j].a_ == acceptors[i]) numHBondsConnToAcceptors++;
            }
        }
        numHBondsConnToAcceptorsArray.emplace_back(numHBondsConnToAcceptors);

        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Print results
    printf("\r\n");
    for(int i=0; i<HBondCriteria.size(); i++)
    {
        fprintf(stdOut_, "\t[%i] -> D:%s, H:%s, A:%s, d-h:%g, h-a:%g, d-a:%g, d-h...a:%g\r\n", i+1,
                HBondCriteria[i].getDonor().data(),
                HBondCriteria[i].getHydrogen().data(),
                HBondCriteria[i].getAcceptor().data(),
                HBondCriteria[i].R_dh(),
                HBondCriteria[i].R_ha(),
                HBondCriteria[i].R_da(),
                HBondCriteria[i].alpha_dha());
    }

    fprintf(stdOut_, "\t%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s%-15s\r\n",
            "Index", "Tot.count", "Num.Donors", "Num.Accept.", "Num.Conn.Don.", "Num.Conn.Acc.", "Frame", "FirstDon.", "FirstHydr.", "FirstAcc");

    for(int i=0; i<HBondTotCounts.size(); i++)
    {
        fprintf(stdOut_, "\t%-15i%-15i%-15i%-15i%-15i%-15i%-15i%-15i%-15i%-15i\r\n", i, HBondTotCounts[i], numDonors[i], numAcceptors[i],
                numHBondsConnToDonorsArray[i], numHBondsConnToAcceptorsArray[i], i*stride+frameFrom, firstHBondInFrame[i].d_, firstHBondInFrame[i].h_, firstHBondInFrame[i].a_);
    }

    return lastError_;
}
