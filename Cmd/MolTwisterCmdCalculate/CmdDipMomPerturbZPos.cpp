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

#include "CmdDipMomPerturbZPos.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdDipMomPerturbZPos::getCmd()
{
    return "dipmomperturbzpos";
}

std::vector<std::string> CCmdDipMomPerturbZPos::getCmdLineKeywords()
{
    return { "dipmomperturbzpos", "chargedmol", "__fromloadedframes__" };
}

std::vector<std::string> CCmdDipMomPerturbZPos::getCmdHelpLines()
{
    return {
                "dipmomperturbzpos <DCD filename> <frame index> <list of atomic IDs> <atoms to displace> <delta Z> <num perturbations> [chargedmol]"
           };
}

std::string CCmdDipMomPerturbZPos::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the dipole moment of the selected molecules <list of atomic IDs> (comma\r\n";
    text+= "\tseparated, no space). The dipolemoment is averaged based on all the defined molecules\r\n";
    text+= "\tin frame index given by <frame index>. The DCD file, <DCD filename>, is used as input,\r\n";
    text+= "\tor the loaded frames if <DCD filename> = __fromloadedframes__. The dipole moment expression\r\n";
    text+= "\tfor neutral molecules are used as default. By sepcifying 'chargedmol', the dipole moment\r\n";
    text+= "\texpression for charged molecules is employed. The dipole moment is perturbed through\r\n";
    text+= "\t<num perturbations> steps, where the <atomd to displace> (comma separated, no space)\r\n";
    text+= "\tare all displaced by the z-position <delta Z>. Note that <atoms to displace> is a lists of\r\n";
    text+= "\tatomic IDs, such as H, O and C7.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Step Px Py Pz Rel.z.pos\r\n";
    text+= "\t2. <step> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Z pos displacement>\r\n";
    text+= "\t3. <step> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Z pos displacement>\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\tN+1. <step> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Z pos displacement>\r\n";
    text+= "\twhere N is the number of perturbations.";

    return text;
}

std::string CCmdDipMomPerturbZPos::execute(std::vector<std::string> arguments)
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

    // Get atom types to displace
    std::vector<std::string> atomsToMove;
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsToMove = CASCIIUtility::getWords(text, ",");

    // Get the z-distance to displace the atoms
    text = CASCIIUtility::getArg(arguments, arg++);
    double displacement = atof(text.data());

    // Get the number of steps
    text = CASCIIUtility::getArg(arguments, arg++);
    int numberOfSteps = atoi(text.data());

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

    // List all atoms to move
    std::vector<int> indicesAtomsToMove;
    for(int i=0; i<(int)atomsToMove.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsToMove[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesAtomsToMove.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // Retrieve DCD at selected frame
    if(useDCD) dcdFile.gotoRecord(frame);

    // Print header
    fprintf(stdOut_, "%-20s%-20s%-20s%-20s%-20s\r\n", "Step", "Px", "Py", "Pz", "Rel.z.pos");

    // Move all atoms that were selected for displacement, stepwize. At each step, print the dipole moment.
    CDCDTools dcdTools(state_, stdOut_);
    double dz = displacement / double(numberOfSteps);
    for(int i=0; i<numberOfSteps; i++)
    {
        // Measure the dipole moment and print it
        C3DVector P, Rc;
        if(useDCD) P = dcdTools.getMoleculeDipoleMoment(indicesDipmomBase, &dcdFile, Rc, chargeNeutralFormulation);
        else P = CMolTwisterStateTools(state_, stdOut_).getMoleculeDipoleMoment(indicesDipmomBase, frame, Rc, chargeNeutralFormulation);
        fprintf(stdOut_, "%-20i%-20.8f%-20.8f%-20.8f%-20.8f\r\n", i, P.x_, P.y_, P.z_, static_cast<double>(i) * dz);

        for(int j=0; j<static_cast<int>(indicesAtomsToMove.size()); j++)
        {
            // Get position of atom
            C3DVector r;
            if(useDCD) r = dcdFile.getCoordinate(indicesAtomsToMove[j]);
            else r = state_->atoms_[indicesAtomsToMove[j]]->r_[frame];

            // Move the atom one step
            r.z_+= dz;
            if(useDCD)
            {
                dcdFile.setCoordinate(indicesAtomsToMove[j], r);
            }
            else
            {
                state_->atoms_[indicesAtomsToMove[j]]->r_[frame] = r;
            }
        }
    }

    return lastError_;
}
