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

#include "CmdEnergyOfTranslation.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../MDFF/MolTwisterMDFFCoulomb.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdEnergyOfTranslation::getCmd()
{
    return "energyoftranslation";
}

std::vector<std::string> CCmdEnergyOfTranslation::getCmdLineKeywords()
{
    return { "energyoftranslation" };
}

std::vector<std::string> CCmdEnergyOfTranslation::getCmdHelpLines()
{
    return {
                "energyoftranslation <xs> <ys> <zs> <xe> <ye> <ze> <steps>"
           };
}

std::string CCmdEnergyOfTranslation::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the total non-bonded energy of the system shown in the current frame for repeated\r\n";
    text+= "\ttranslations of the selected atoms in the current frame. The selected atoms are translated along\r\n";
    text+= "\tfrom start vector, (<xs>, <ys>, <zs>) to the end vector, (<xe> <ye> <ze>), in <step> steps. The\r\n";
    text+= "\tnon-bonded force fields must be loaded or created before running this command.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Position ENon-bond\r\n";
    text+= "\t2. <x> <y> <z> <total system energy (kJ/mol)>\r\n";
    text+= "\t3. <x> <y> <z> <total system energy (kJ/mol)>\r\n";
    text+= "\t              .\r\n";
    text+= "\t              .\r\n";
    text+= "\t              .\r\n";
    text+= "\tN+1. <x> <y> <z> <total system energy (kJ/mol)>\r\n";
    text+= "\twhere N is the number of steps. The position (<x>, <y>, <z>) are the positions between the\r\n";
    text+= "\tstart vector and the end vector.";

    return text;
}

std::string CCmdEnergyOfTranslation::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<int> atomIndicesToInclude;
    C3DVector vecStart, vecEnd;
    std::string text;
    int steps;


    // Retrieve start point and end point of translation, as well as number of steps
    text = CASCIIUtility::getArg(arguments, arg++);
    vecStart.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecStart.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecStart.z_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecEnd.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecEnd.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecEnd.z_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    steps = atoi(text.data());
    if(steps < 2)
    {
        lastError_ = "Error: number of steps must be greater than one!";
        return lastError_;
    }


    // Retrieve atoms to include in move
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i]->isSelected())
            atomIndicesToInclude.emplace_back(i);
    }


    // Find center of mass for atoms in move and calculate translation vector
    // onto new positions for the atoms
    C3DVector vecCOM = CMolTwisterStateTools(state_, stdOut_).getCenterOfMass(atomIndicesToInclude, state_->currentFrame_);
    C3DVector vecTrans = vecStart - vecCOM;


    // Find starting positions of all atoms included in move,
    // and also store old postions
    std::vector<C3DVector> translatedPositions;
    std::vector<C3DVector> oldPositions;
    translatedPositions.resize(atomIndicesToInclude.size());
    oldPositions.resize(atomIndicesToInclude.size());
    for(int i=0; i<atomIndicesToInclude.size(); i++)
    {
        CAtom* atom = state_->atoms_[atomIndicesToInclude[i]].get();
        C3DVector vecCurrPos = atom->r_[state_->currentFrame_];
        oldPositions[i] = vecCurrPos;
        translatedPositions[i] = vecCurrPos + vecTrans;
    }


    // Find vector to translate for each step
    vecTrans = (vecEnd - vecStart)*(1.0 / double(steps - 1));


    // Step from vStart to vEnd and estimate total non-bonded energy at each step
    printf("\r\n");
    fprintf(stdOut_, "\t%-30s%-15s\r\n", "Position", "ENon-bond");
    char posString[128];
    CMDFFCoulomb coulomb(state_);
    pb.beginProgress("Calculate total energy as function of position");
    for(int n=0; n<steps; n++)
    {
        // Translate atoms to move
        for(int i=0; i<translatedPositions.size(); i++)
        {
            state_->atoms_[atomIndicesToInclude[i]]->r_[state_->currentFrame_] = translatedPositions[i];
        }


        // Calculate total non-bonded energy
        double ETot = 0.0;
        for(int I=0; I<(state_->atoms_.size()-1); I++)
        {
            C3DVector rI = state_->atoms_[I]->r_[state_->currentFrame_];
            double qI = state_->atoms_[I]->Q_;
            std::string IDI = state_->atoms_[I]->getID();
            for(int J=I+1; J<state_->atoms_.size(); J++)
            {
                C3DVector rJ = state_->atoms_[J]->r_[state_->currentFrame_];
                double qJ = state_->atoms_[J]->Q_;
                std::string IDJ = state_->atoms_[J]->getID();

                ETot+= coulomb.calcPotentialBetween(rI, rJ, qI, qJ);

                std::shared_ptr<std::vector<int>> ffIndices = state_->mdFFNonBondedList_.indexFromNames(IDI, IDJ);
                for(int k=0; k<ffIndices->size(); k++)
                {
                    ETot+= state_->mdFFNonBondedList_.get((*ffIndices)[k])->calcPotential(rI, rJ);
                }
            }
        }


        // Report total non-boned energy
        C3DVector vecCurrPos = vecStart + vecTrans*(double(n));
        sprintf(posString, "%.4f %.4f %.4f", vecCurrPos.x_, vecCurrPos.y_, vecCurrPos.z_);
        fprintf(stdOut_, "\t%-30s%-15.8f\r\n", posString, ETot);


        // Update positions for next step
        for(int i=0; i<translatedPositions.size(); i++)
            translatedPositions[i]+= vecTrans;

        pb.updateProgress(n, steps);
    }
    pb.endProgress();


    // Restore old positons
    for(int i=0; i<oldPositions.size(); i++)
    {
        state_->atoms_[atomIndicesToInclude[i]]->r_[state_->currentFrame_] = oldPositions[i];
    }

    return lastError_;
}
