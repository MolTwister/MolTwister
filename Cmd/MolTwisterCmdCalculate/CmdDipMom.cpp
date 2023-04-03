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

#include "CmdDipMom.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdDipMom::getCmd()
{
    return "dipmom";
}

std::vector<std::string> CCmdDipMom::getCmdLineKeywords()
{
    return { "dipmom", "__fromloadedframes__", "chargedmol" };
}

std::vector<std::string> CCmdDipMom::getCmdHelpLines()
{
    return {
                "dipmom <DCD filename> <frame index> <list of atomic IDs (comma sep., no space)> [chargedmol]"
           };
}

std::string CCmdDipMom::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the dipole moment of the selected molecules <list of atomic IDs>. The\r\n";
    text+= "\tdipolemoment is averaged based on all the defined molecules in frame index given by\r\n";
    text+= "\t<frame index>. Either one can choose a DCD file through <DCD filename> as input, or\r\n";
    text+= "\tit is possible to use the loaded frames as input by letting <DCD filename> = \r\n";
    text+= "\t__fromloadedframes__. The dipole moment expression for neutral molecules are used as\r\n";
    text+= "\tdefault. By sepcifying 'chargedmol', the dipole moment expression for charged molecules\r\n";
    text+= "\tis employed. The output is defined below.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tPx Py Pz\r\n";
    text+= "\t<dipole moment x-component> <dipole moment y-component> <dipole moment z-component>";

    return text;
}

std::string CCmdDipMom::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CDCDFile dcdFile;

    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);

    bool useDCD = false;
    if(text != "__fromloadedframes__")
    {
        if(!dcdFile.open(text))
        {
            lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
            return lastError_;
        }

        useDCD = true;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    int frame = atoi(text.data());
    if(useDCD)
    {
        if((frame < 0) || (frame >= dcdFile.getNumRecords()))
        {
            lastError_ = "Error: Frame is outside the range of available number of frames (note: index should be zero indexed)";
            return lastError_;
        }
    }
    else
    {
        if(state_->atoms_.size() < 1)
        {
            lastError_ = "Error: no atoms found";
            return lastError_;
        }
        if((frame < 0) || (frame >= state_->atoms_[0]->r_.size()))
        {
            lastError_ = "Error: Frame is outside the range of available number of frames (note: index should be zero indexed)";
            return lastError_;
        }
    }

    // Get atom types to calculate the dipole moment for
    std::vector<std::string> atomsDipmomBase;
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsDipmomBase = CASCIIUtility::getWords(text, ",");

    // Are we going to use charge neutral formulation: sum_i q_i r_i or
    // are we using charged molecules: sum_i q_i (r_i - R_c)?
    bool chargeNeutralFormulation = true;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "chargedmol")
        chargeNeutralFormulation = false;
    else arg--;

    // List all selected atom indices to calculate the dipole moment for
    std::vector<int> indicesDipmomBase;
    for(int i=0; i<(int)atomsDipmomBase.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsDipmomBase[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesDipmomBase.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // Retrieve DCD at selected frame
    if(useDCD) dcdFile.gotoRecord(frame);

    // Print header
    fprintf(stdOut_, "%-20s%-20s%-20s\r\n", "Px", "Py", "Pz");

    // Measure the dipole moment and print it
    C3DVector P, Rc;
    if(useDCD) P = CDCDTools(state_, stdOut_).getMoleculeDipoleMoment(indicesDipmomBase, &dcdFile, Rc, chargeNeutralFormulation);
    else P = CMolTwisterStateTools(state_, stdOut_).getMoleculeDipoleMoment(indicesDipmomBase, frame, Rc, chargeNeutralFormulation);
    fprintf(stdOut_, "%-20.8f%-20.8f%-20.8f\r\n", P.x_, P.y_, P.z_);

    return lastError_;
}
