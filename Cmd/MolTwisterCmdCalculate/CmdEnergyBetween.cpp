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

#include "CmdEnergyBetween.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../MDFF/MolTwisterMDFFCoulomb.h"

std::string CCmdEnergyBetween::getCmd()
{
    return "energybetween";
}

std::vector<std::string> CCmdEnergyBetween::getCmdLineKeywords()
{
    return { "energybetween", "nonbonded", "coulomb", "bond", "angle", "dihedral" };
}

std::vector<std::string> CCmdEnergyBetween::getCmdHelpLines()
{
    return {
                "energybetween nonbonded <comma sep. list of atom IDs> <FF index (zero based)>",
                "energybetween coulomb <comma sep. list of atom IDs>",
                "energybetween bond <comma sep. list of atom IDs> <FF index (zero based)>",
                "energybetween angle <comma sep. list of atom IDs> <FF index (zero based)>",
                "energybetween dihedral <comma sep. list of atom IDs> <FF index (zero based)>"
           };
}

std::string CCmdEnergyBetween::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the nonbonded, Coulomb or bond energy (in kJ/mol) betwee two atoms, based on the\r\n";
    text+= "\tconfigured force fields, or the angular energy between three atoms or the dihedral energy\r\n";
    text+= "\tbetween four atoms. The energies are calculated based on the atoms shown in the current frame.\r\n";
    text+= "\tThe indices of the various force fields that have been created or loaded is available through\r\n";
    text+= "\tthe 'list' command. An index into these lists, for the appropriate force field type, is specified\r\n";
    text+= "\tthrough the <FF index (zero based)> parameter. The energy is calculated between the atoms given\r\n";
    text+= "\tin <comma sep. list of atom IDs> (e.g., H, O, C7, where the list must be enetered withour spece).";

    return text;
}

std::string CCmdEnergyBetween::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text, energyType;
    std::vector<std::string> atomsToInclude;


    // Get force type
    energyType = CASCIIUtility::getArg(arguments, arg++);

    // Find atom indices to calculate energy between
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsToInclude = CASCIIUtility::getWords(text, ",");

    // Calculate forces
    if(energyType == "nonbonded")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for non-bonded energy calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);
        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFNonBondedList_.size()))
        {
            lastError_ = std::string("Error: non-bonded force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        if(!state_->mdFFNonBondedList_.get(ffIndex)->isPotentialCalcAvailable())
        {
            lastError_ = std::string("Error: non-bonded force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of energy calculations!");
            return lastError_;
        }

        double U = state_->mdFFNonBondedList_.get(ffIndex)->calcPotential(r1, r2);

        printf("\r\n");
        fprintf(stdOut_, "\tEnergy (%s) between %s and %s (index %i and %i) at pos (%g,%g,%g) and (%g,%g,%g): U = %g\r\n", state_->mdFFNonBondedList_.get(ffIndex)->getFFType().data(), atomName1.data(), atomName2.data(), index1, index2, r1.x_, r1.y_, r1.z_, r2.x_, r2.y_, r2.z_, U);
    }
    else if(energyType == "coulomb")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for Coulomb energy calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        double q1 = state_->atoms_[index1]->Q_;
        double q2 = state_->atoms_[index2]->Q_;
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        CMDFFCoulomb Coulomb(state_);
        double U = Coulomb.calcPotentialBetween(r1, r2, q1, q2);

        printf("\r\n");
        fprintf(stdOut_, "\tEnergy (Coulomb) between %s and %s (index %i and %i) at pos (%g,%g,%g) and (%g,%g,%g): U = %g\r\n", atomName1.data(), atomName2.data(), index1, index2, r1.x_, r1.y_, r1.z_, r2.x_, r2.y_, r2.z_, U);
    }
    else if(energyType == "bond")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for bond energy calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFBondList_.size()))
        {
            lastError_ = std::string("Error: bond force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        if(!state_->mdFFBondList_.get(ffIndex)->isPotentialCalcAvailable())
        {
            lastError_ = std::string("Error: bond force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of energy calculations!");
            return lastError_;
        }

        double U = state_->mdFFBondList_.get(ffIndex)->calcPotential(r1, r2);

        printf("\r\n");
        fprintf(stdOut_, "\tEnergy (%s) between %s and %s (index %i and %i) at pos (%g,%g,%g) and (%g,%g,%g): U = %g\r\n", state_->mdFFBondList_.get(ffIndex)->getFFType().data(), atomName1.data(), atomName2.data(), index1, index2, r1.x_, r1.y_, r1.z_, r2.x_, r2.y_, r2.z_, U);
    }
    else if(energyType == "angle")
    {
        if(atomsToInclude.size() != 3)
        {
            lastError_ = "Error: only three atoms can be specified for angle energy calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        int index3 = atoi(atomsToInclude[2].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }
        if((index3 < 0) || (index3 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index3) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFAngleList_.size()))
        {
            lastError_ = std::string("Error: angle force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        C3DVector r3 = state_->atoms_[index3]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();
        std::string atomName3 = state_->atoms_[index3]->getID();

        if(!state_->mdFFAngleList_.get(ffIndex)->isPotentialCalcAvailable())
        {
            lastError_ = std::string("Error: angle force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of energy calculations!");
            return lastError_;
        }

        double U = state_->mdFFAngleList_.get(ffIndex)->calcPotential(r1, r2, r3);

        printf("\r\n");
        fprintf(stdOut_, "\tEnergy (%s) between %s, %s and %s (index %i, %i and %i) at pos (%g,%g,%g), (%g,%g,%g) and (%g,%g,%g): U = %g\r\n", state_->mdFFAngleList_.get(ffIndex)->getFFType().data(), atomName1.data(), atomName2.data(), atomName3.data(), index1, index2, index3, r1.x_, r1.y_, r1.z_, r2.x_, r2.y_, r2.z_, r3.x_, r3.y_, r3.z_, U);
    }
    else if(energyType == "dihedral")
    {
        if(atomsToInclude.size() != 4)
        {
            lastError_ = "Error: only four atoms can be specified for dihedral energy calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        int index3 = atoi(atomsToInclude[2].data());
        int index4 = atoi(atomsToInclude[3].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }
        if((index3 < 0) || (index3 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index3) + std::string(" does not exist!");
            return lastError_;
        }
        if((index4 < 0) || (index4 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index4) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFDihList_.size()))
        {
            lastError_ = std::string("Error: dihedral force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        C3DVector r3 = state_->atoms_[index3]->r_[state_->currentFrame_];
        C3DVector r4 = state_->atoms_[index4]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();
        std::string atomName3 = state_->atoms_[index3]->getID();
        std::string atomName4 = state_->atoms_[index4]->getID();

        if(!state_->mdFFDihList_.get(ffIndex)->isPotentialCalcAvailable())
        {
            lastError_ = std::string("Error: dihedral force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of energy calculations!");
            return lastError_;
        }

        double U = state_->mdFFDihList_.get(ffIndex)->calcPotential(r1, r2, r3, r4);

        printf("\r\n");
        fprintf(stdOut_, "\tEnergy (%s) between %s, %s, %s and %s (index %i, %i, %i and %i) at pos (%g,%g,%g), (%g,%g,%g), (%g,%g,%g) and (%g,%g,%g): U = %g\r\n", state_->mdFFDihList_.get(ffIndex)->getFFType().data(), atomName1.data(), atomName2.data(), atomName3.data(), atomName4.data(), index1, index2, index3, index4, r1.x_, r1.y_, r1.z_, r2.x_, r2.y_, r2.z_, r3.x_, r3.y_, r3.z_, r4.x_, r4.y_, r4.z_, U);
    }
    else
    {
        lastError_ = "Error: expected 'nonbonded', 'bonded', 'angle' or 'dihedral'!";
        return lastError_;
    }

    return lastError_;
}
