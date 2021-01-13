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

#include "CmdDipMomPerturbCharge.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdDipMomPerturbCharge::getCmd()
{
    return "dipmomperturbcharge";
}

std::vector<std::string> CCmdDipMomPerturbCharge::getCmdLineKeywords()
{
    return { "dipmomperturbcharge", "chargedmol" };
}

std::vector<std::string> CCmdDipMomPerturbCharge::getCmdHelpLines()
{
    return {
                "dipmomperturbcharge <DCD filename> <frame index> <list of atomic IDs> <positive charges> <negative charges> <num perturbations> <delta Q> [chargedmol]"
           };
}

std::string CCmdDipMomPerturbCharge::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the dipole moment of the selected molecules <list of atomic IDs> (comma\r\n";
    text+= "\tseparated, no space). The dipolemoment is averaged based on all the defined molecules\r\n";
    text+= "\tin frame index given by <frame index>. The DCD file, <DCD filename>, is used as input.\r\n";
    text+= "\tThe dipole moment expression for neutral molecules are used as default. By sepcifying\r\n";
    text+= "\t'chargedmol', the dipole moment expression for charged molecules is employed. The dipole\r\n";
    text+= "\tmoment is perturbed through <num perturbations> steps, where the positive and negative\r\n";
    text+= "\tcharges are all increased and decreased, respectively, by <delta Q> spread over the\r\n";
    text+= "\trespective groups of charges, thus preserving charge neutrality. The greatest perturbation\r\n";
    text+= "\tis calculated first. Note that the lists of positive and negative charges are lists of\r\n";
    text+= "\tatomic IDs, such as H, O and C7 (comma separated, no space).\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. QPerturb+ QPerturb- Px Py Pz\r\n";
    text+= "\t2. <positive perturbation> <negative perturbation> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component>\r\n";
    text+= "\t3. <positive perturbation> <negative perturbation> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component>\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\tN+1. <positive perturbation> <negative perturbation> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component>\r\n";
    text+= "\twhere N is the number of perturbations.";

    return text;
}

std::string CCmdDipMomPerturbCharge::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CDCDFile dcdFile;

    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    int frame = atoi(text.data());
    if((frame < 0) || (frame >= dcdFile.getNumRecords()))
    {
        lastError_ = "Error: Frame is outside the range of available number of frames (note: index should be zero indexed)";
        return lastError_;
    }

    // Get atom types to calculate the dipole moment for
    std::vector<std::string> atomsDipmomBase;
    text = CASCIIUtility::getArg(arguments, arg++);

    CASCIIUtility::removeWhiteSpace(text);
    atomsDipmomBase = CASCIIUtility::getWords(text, ",");

    // Get atom types with positive charge
    std::vector<std::string> atomsPositive;
    text = CASCIIUtility::getArg(arguments, arg++);

    CASCIIUtility::removeWhiteSpace(text);
    atomsPositive = CASCIIUtility::getWords(text, ",");

    // Get atom types with negative charge
    std::vector<std::string> atomsNegative;
    text = CASCIIUtility::getArg(arguments, arg++);

    CASCIIUtility::removeWhiteSpace(text);
    atomsNegative = CASCIIUtility::getWords(text, ",");

    // Get number of steps to perturb the charge, as well as the step size
    text = CASCIIUtility::getArg(arguments, arg++);
    int numPerturbations = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double deltaQ = atof(text.data());

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

    // List all selected atom indices associated with a positive charge
    std::vector<int> indicesPositive;
    for(int i=0; i<(int)atomsPositive.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsPositive[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesPositive.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // List all selected atom indices associated with a negative charge
    std::vector<int> indicesNegative;
    for(int i=0; i<(int)atomsNegative.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsNegative[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesNegative.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // Retrieve DCD at selected frame
    dcdFile.gotoRecord(frame);

    // Print header
    fprintf(stdOut_, "%-20s%-20s%-20s%-20s%-20s\r\n", "QPerturb+", "QPerturb-", "Px", "Py", "Pz");

    // Calculate dipole moment while incrementally perturbing the selected positive
    // and negative charges, such that the system remains charge neutral. First increase
    // positive charges and decrease negative ones (positive Qperturb). Then move further
    // and increase the negative charges and decrease the positive ones (negative Qperturb).
    C3DVector P, Rc;
    CDCDTools dcdTools(state_, stdOut_);
    double deltaQPositive = deltaQ / double(indicesPositive.size());
    double deltaQNegative = deltaQ / double(indicesNegative.size());
    double qPerturbPositive = deltaQPositive * double(numPerturbations);
    double qPerturbNegative = deltaQNegative * double(numPerturbations);
    for(int i=0; i<(2*numPerturbations); i++)
    {
        // Perturb selected atoms with positive charge
        for(int j=0; j<indicesPositive.size(); j++)
            state_->atoms_[indicesPositive[j]]->Q_+= qPerturbPositive;

        // Perturb selected atoms with negative charge
        for(int j=0; j<indicesNegative.size(); j++)
            state_->atoms_[indicesNegative[j]]->Q_-= qPerturbNegative;

        // Measure the dipole moment and print it
        P = dcdTools.getMoleculeDipoleMoment(indicesDipmomBase, &dcdFile, Rc, chargeNeutralFormulation);
        fprintf(stdOut_, "%-20.8f%-20.8f%-20.8f%-20.8f%-20.8f\r\n", qPerturbPositive, -qPerturbNegative, P.x_, P.y_, P.z_);

        // Restore the original charges
        for(int j=0; j<indicesPositive.size(); j++)
            state_->atoms_[indicesPositive[j]]->Q_-= qPerturbPositive;
        for(int j=0; j<indicesNegative.size(); j++)
            state_->atoms_[indicesNegative[j]]->Q_+= qPerturbNegative;

        qPerturbPositive-= deltaQPositive;
        qPerturbNegative-= deltaQNegative;
    }

    return lastError_;
}
