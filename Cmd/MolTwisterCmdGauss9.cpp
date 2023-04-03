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

#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "MolTwisterCmdMeasure.h"
#include "MolTwisterCmdGauss9.h"
#include "Tools/MolecularTools.h"

void CCmdGauss9::onAddKeywords()
{
    addKeyword("gauss9");
    addKeyword("dihedralrot");
    addKeyword("rot");
    addKeyword("qmspec");
    addKeyword("genxyzfrominput");
    addKeyword("relaxed");
    addKeyword("genoptefromoutput");
    addKeyword("genoptxyzfromoutput");
    addKeyword("genxyzfromoutput");
    addKeyword("anglerot");
    addKeyword("bondstretch");
    addKeyword("stretch");
}

std::string CCmdGauss9::getHelpString() const
{
    std::string text;

    text+= "\tUsage: gauss9 dihedralrot id <atom index 1> <atom index 2> <atom index 3> <atom index41> rot <start angle> <end angle> <angular step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 dihedralrot var <variable name> rot <start angle> <end angle> <angular step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 anglerot id <atom index 1> <atom index 2> <atom index 3> rot <start angle> <end angle> <angular step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 anglerot var <variable name> rot <start angle> <end angle> <angular step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 bondstretch id <atom index 1> <atom index 2> stretch <start dist> <end dist> <dist step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 bondstretch var <variable name> stretch <start dist> <end dist> <dist step size> qmspec <base set> <charge> <spin multiplicity> [relaxed]\r\n";
    text+= "\t       gauss9 genxyzfrominput <gaussian input file name>\r\n";
    text+= "\t       gauss9 genoptefromoutput <gaussian output file name>\r\n";
    text+= "\t       gauss9 genoptxyzfromoutput <gaussian output file name>\r\n";
    text+= "\t       gauss9 genxyzfromoutput <gaussian output file name>\r\n";
    text+= "\r\n";
    text+= "\tThis command can execute several operations that are useful for the Gaussian9 software package. These are listed in the following:\r\n";
    text+= "\t* dihedralrot: will generate a Gaussian9 input script that executes electronic structure calculations on the loaded system\r\n";
    text+= "\t               for each rotation of a dihedral, either specified through 4 atomic indices or through the same indices stored\r\n";
    text+= "\t               in a variable.\r\n";
    text+= "\t* anglerot:    will generate a Gaussian9 input script that executes electronic structure calculations on the loaded system\r\n";
    text+= "\t               for each rotation of an angle, either specified through 3 atomic indices or through the same indices stored\r\n";
    text+= "\t               in a variable.\r\n";
    text+= "\t* bondstretch: will generate a Gaussian9 input script that executes electronic structure calculations on the loaded system\r\n";
    text+= "\t               for each atomic distance separation, either specified through 2 atomic indices or through the same indices\r\n";
    text+= "\t               stored in a variable.\r\n";
    text+= "\t* genxyzfrominput: Extracts the atom names and coordinates from a Gaussian9 input file and outputs a XYZ-file formatted text.\r\n";
    text+= "\t* genoptefromoutput: Extracts the optimization energies from a Gaussian9 output file and produces a list of these.\r\n";
    text+= "\t* genoptxyzfromoutput: Extracts the optimized structures from a Gaussian9 output file and produces XYZ-formatted text outputs.\r\n";
    text+= "\t* genxyzfromoutput: Extracts the input structure from a Gaussian9 output file and produces XYZ-formatted text outputs.";

    return text;
}

void CCmdGauss9::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
     
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "dihedralrot")
    {
        parseDihedralrotCommand(commandLine, arg);
    }

    else if(text == "anglerot")
    {
        parseAnglerotCommand(commandLine, arg);
    }

    else if(text == "bondstretch")
    {
        parseBondstretchCommand(commandLine, arg);
    }
    
    else if(text == "genxyzfrominput")
    {
        parseGenxyzfrominputCommand(commandLine, arg);
    }
    
    else if(text == "genoptefromoutput")
    {
        parseGenoptefromoutputCommand(commandLine, arg);
    }

    else if(text == "genoptxyzfromoutput")
    {
        parseGenoptxyzfromoutputCommand(commandLine, arg);
    }

    else if(text == "genxyzfromoutput")
    {
        parseGenxyzfromoutputCommand(commandLine, arg);
    }
    
    else
    {
        printf("Syntax Error: Second argument should specify the kind of Gaussian script to generate!");
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1) state_->view3D_->requestUpdate(true);
        else state_->view3D_->requestUpdate(false);
    }
}

void CCmdGauss9::parseDihedralrotCommand(std::string commandLine, int& arg)
{
    std::string text;
    CAtom* atom1Ptr = nullptr;
    CAtom* atom2Ptr = nullptr;
    CAtom* atom3Ptr = nullptr;
    CAtom* atom4Ptr = nullptr;
    bool foundIndex = false;
    int atomIndex1=0, atomIndex2=0, atomIndex3=0, atomIndex4=0;
      
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex1 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex2 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex3 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex4 = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeDihedral))
        {
            CVarDihedral* p = (CVarDihedral*)varPtr;
            
            atomIndex1 = p->atomIndex1_;
            atomIndex2 = p->atomIndex2_;
            atomIndex3 = p->atomIndex3_;
            atomIndex4 = p->atomIndex4_;
            foundIndex = true;
        }
        else
        {
            printf("Variable missing or not of type 'dihedral'!");
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'id' or 'var'!");
    }
    
    if(foundIndex)
    {
        if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size())
           && (atomIndex3 < state_->atoms_.size()) && (atomIndex4 < state_->atoms_.size()))
        {
            atom1Ptr = state_->atoms_[atomIndex1].get();
            atom2Ptr = state_->atoms_[atomIndex2].get();
            atom3Ptr = state_->atoms_[atomIndex3].get();
            atom4Ptr = state_->atoms_[atomIndex4].get();
            if(atom1Ptr && atom2Ptr && atom3Ptr && atom4Ptr)
            {
                text = CASCIIUtility::getWord(commandLine, arg++);
                if(text == "rot")
                {
                    double startAngle, endAngle, step;

                    text = CASCIIUtility::getWord(commandLine, arg++);
                    startAngle = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    endAngle = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    step = atof(text.data());
                    if(step == 0.0) step = (endAngle - startAngle) / 10.0;
                    if(startAngle > endAngle) endAngle = startAngle + 360.0;

                    text = CASCIIUtility::getWord(commandLine, arg++);
                    if(text == "qmspec")
                    {
                        std::string baseSetSpec;
                        double charge;
                        bool relaxed = false;
                        int spinMul, index = 0;

                        baseSetSpec = CASCIIUtility::getWord(commandLine, arg++);
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        charge = atof(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        spinMul = atoi(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        if(text == "relaxed") relaxed = true;

                        if(relaxed)
                        {
                            // Generate gaussian script that runs an optimization at each dihedral
                            // angle, but keeps the rotating dihedral fixed during optimization
                            // (i.e. relaxed scan)
                            std::string str;
                            str.resize(512);
                            sprintf((char*)str.data(), "Relaxed dihedral scan in range Phi=[%g, %g]deg, Step=%g", startAngle, endAngle, step);
                            printGaussianModRedScript(baseSetSpec.data(), charge, spinMul, str.data(),
                                                      atomIndex1, atomIndex2, atomIndex3, atomIndex4, startAngle, endAngle, step);
                        }
                        else
                        {
                            // Generate gaussian script that rotates the dihedral as well as all
                            // atoms bonded to dihedral and perform SP (single point) energy calc. 
                            // (i.e. rigid scan)
                            if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                            {
                                for(double angle=startAngle; angle<=endAngle; angle+=step)
                                {
                                    CMolecularTools::modDihedralTo(atom1Ptr, atom2Ptr, atom3Ptr, atom4Ptr, angle * M_PI/180.0, state_->getCurrFrameIndex());
                                    
                                    if(index > 0) fprintf(stdOut_, "\t--Link1--\r\n\r\n");
                                    index++;
                                    std::string str;
                                    str.resize(512);
                                    sprintf((char*)str.data(), "Single point energy calculation Phi=%g deg", angle);
                                    printGaussianSPScript(baseSetSpec.data(), charge, spinMul, str.data());
                                }
                                state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                            }
                        }
                    }
                    else
                    {
                        printf("Syntax: Error: Seventh argument should be 'qmspec', followed by <gaussian basis spec.> <charge> <spin mul.>!");
                    }
                }
                else
                {
                    printf("Syntax: Error: Fifth argument should be 'rot', followed by <start> <end> <step>!");
                }
            }
            else
            {
                printf("Could not find requested atoms!");
            }
        }
        else
        {
            printf("Invalid atom indices!");
        }
    }
}

void CCmdGauss9::parseAnglerotCommand(std::string commandLine, int& arg)
{
    std::string text;
    CAtom* atom1Ptr = nullptr;
    CAtom* atom2Ptr = nullptr;
    CAtom* atom3Ptr = nullptr;
    bool foundIndex = false;
    int atomIndex1=0, atomIndex2=0, atomIndex3=0;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex1 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex2 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex3 = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeAngle))
        {
            CVarAngle* p = (CVarAngle*)varPtr;
            
            atomIndex1 = p->atomIndex1_;
            atomIndex2 = p->atomIndex2_;
            atomIndex3 = p->atomIndex3_;
            foundIndex = true;
        }
        else
        {
            printf("Variable missing or not of type 'angle'!");
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'id' or 'var'!");
    }
    
    if(foundIndex)
    {
        if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size())
           && (atomIndex3 < state_->atoms_.size()))
        {
            atom1Ptr = state_->atoms_[atomIndex1].get();
            atom2Ptr = state_->atoms_[atomIndex2].get();
            atom3Ptr = state_->atoms_[atomIndex3].get();
            if(atom1Ptr && atom2Ptr && atom3Ptr)
            {
                text = CASCIIUtility::getWord(commandLine, arg++);
                if(text == "rot")
                {
                    double startAngle, endAngle, step;
                    
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    startAngle = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    endAngle = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    step = atof(text.data());
                    if(step == 0.0) step = (endAngle - startAngle) / 10.0;
                    if(startAngle > endAngle) endAngle = startAngle + 360.0;
                    
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    if(text == "qmspec")
                    {
                        std::string baseSetSpec;
                        double charge;
                        bool relaxed = false;
                        int spinMul, index = 0;
                        
                        baseSetSpec = CASCIIUtility::getWord(commandLine, arg++);
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        charge = atof(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        spinMul = atoi(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        if(text == "relaxed") relaxed = true;
                        
                        if(relaxed)
                        {
                            // Generate gaussian script that runs an optimization at each covalent
                            // angle, but keeps the rotating angle fixed during optimization
                            // (i.e. relaxed scan)
                            std::string str;
                            str.resize(512);
                            sprintf((char*)str.data(), "Relaxed covalent angle scan in range Theta=[%g, %g]deg, Step=%g", startAngle, endAngle, step);
                            printGaussianModRedScript(baseSetSpec.data(), charge, spinMul, str.data(),
                                                      atomIndex1, atomIndex2, atomIndex3, startAngle, endAngle, step);
                        }
                        else
                        {
                            // Generate gaussian script that rotates the covalent angle as well as all
                            // atoms bonded to angle and perform SP (single point) energy calc. 
                            // (i.e. rigid scan)
                            if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                            {
                                for(double dAngle=startAngle; dAngle<=endAngle; dAngle+=step)
                                {
                                    CMolecularTools::modAngleTo(atom1Ptr, atom2Ptr, atom3Ptr, dAngle * M_PI/180.0, state_->getCurrFrameIndex());
                                    
                                    if(index > 0) fprintf(stdOut_, "\t--Link1--\r\n\r\n");
                                    index++;
                                    std::string str;
                                    str.resize(512);
                                    sprintf((char*)str.data(), "Single point energy calculation Theta=%g deg", dAngle);
                                    printGaussianSPScript(baseSetSpec.data(), charge, spinMul, str.data());
                                }
                                state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                            }
                        }
                    }
                    else
                    {
                        printf("Syntax: Error: Seventh argument should be 'qmspec', followed by <gaussian basis spec.> <charge> <spin mul.>!");
                    }
                }
                else
                {
                    printf("Syntax: Error: Fifth argument should be 'rot', followed by <start> <end> <step>!");
                }
            }
            else
            {
                printf("Could not find requested atoms!");
            }
        }
        else
        {
            printf("Invalid atom indices!");
        }
    }
}

void CCmdGauss9::parseBondstretchCommand(std::string commandLine, int& arg)
{
    std::string text;
    CAtom* atom1Ptr = nullptr;
    CAtom* atom2Ptr = nullptr;
    bool foundIndex = false;
    int atomIndex1=0, atomIndex2=0;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex1 = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex2 = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeBond))
        {
            CVarBond* p = (CVarBond*)varPtr;
            
            atomIndex1 = p->atomIndex1_;
            atomIndex2 = p->atomIndex2_;
            foundIndex = true;
        }
        else
        {
            printf("Variable missing or not of type 'angle'!");
        }
    }
    else
    {
        printf("Syntax Error: Third argument should be 'id' or 'var'!");
    }
    
    if(foundIndex)
    {
        if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size()))
        {
            atom1Ptr = state_->atoms_[atomIndex1].get();
            atom2Ptr = state_->atoms_[atomIndex2].get();
            if(atom1Ptr && atom2Ptr)
            {
                text = CASCIIUtility::getWord(commandLine, arg++);
                if(text == "stretch")
                {
                    double startDist, endDist, step;
                    
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    startDist = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    endDist = atof(text.data());
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    step = atof(text.data());
                    if(step == 0.0) step = (endDist - startDist) / 10.0;
                    
                    text = CASCIIUtility::getWord(commandLine, arg++);
                    if(text == "qmspec")
                    {
                        std::string baseSetSpec;
                        double charge;
                        bool relaxed = false;
                        int spinMul, index = 0;
                        
                        baseSetSpec = CASCIIUtility::getWord(commandLine, arg++);
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        charge = atof(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        spinMul = atoi(text.data());
                        text = CASCIIUtility::getWord(commandLine, arg++);
                        if(text == "relaxed") relaxed = true;
                        
                        if(relaxed)
                        {
                            // Generate gaussian script that runs an optimization at each bond
                            // length, but keeps the bond length fixed during optimization
                            // (i.e. relaxed scan)
                            std::string str;
                            str.resize(512);
                            sprintf((char*)str.data(), "Relaxed bond length scan in range R=[%g, %g]AA, Step=%g", startDist, endDist, step);
                            printGaussianModRedScript(baseSetSpec.data(), charge, spinMul, str.data(),
                                                      atomIndex1, atomIndex2, startDist, endDist, step);
                        }
                        else
                        {
                            // Generate gaussian script that rotates the covalent angle as well as all
                            // atoms bonded to angle and perform SP (single point) energy calc. 
                            // (i.e. rigid scan)
                            if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                            {
                                for(double dDist=startDist; dDist<=endDist; dDist+=step)
                                {
                                    CMolecularTools::modBondLengthTo(atom1Ptr, atom2Ptr, dDist, state_->getCurrFrameIndex());
                                    
                                    if(index > 0) fprintf(stdOut_, "\t--Link1--\r\n\r\n");
                                    index++;
                                    std::string str;
                                    str.resize(512);
                                    sprintf((char*)str.data(), "Single point energy calculation R=%g AA", dDist);
                                    printGaussianSPScript(baseSetSpec.data(), charge, spinMul, str.data());
                                }
                                state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                            }
                        }
                    }
                    else
                    {
                        printf("Syntax: Error: Seventh argument should be 'qmspec', followed by <gaussian basis spec.> <charge> <spin mul.>!");
                    }
                }
                else
                {
                    printf("Syntax: Error: Fifth argument should be 'stretch', followed by <start> <end> <step>!");
                }
            }
            else
            {
                printf("Could not find requested atoms!");
            }
        }
        else
        {
            printf("Invalid atom indices!");
        }
    }
}

void CCmdGauss9::printGaussianSPScript(std::string baseSetSpec, double charge, int spinMul, std::string heading) const
{
    C3DVector r;
    CAtom* atomPtr;
    
    fprintf(stdOut_, "\t%%Chk = CheckPnt\r\n");
    fprintf(stdOut_, "\t# %s\r\n\r\n\t%s\r\n\r\n\t%g %i\r\n", baseSetSpec.data(), heading.data(), charge, spinMul);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        if(atomPtr)
        {
            std::string ID = atomPtr->getID();
            if(state_->getCurrFrameIndex() < atomPtr->r_.size())
            {
                r = atomPtr->r_[state_->getCurrFrameIndex()];
                fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), r.x_, r.y_, r.z_);
            }
            else
            {
                fprintf(stdOut_, "\t%-5s%-13s%-13s%-13s\r\n", ID.data(), "na", "na", "na");
            }
        }
    }
    fprintf(stdOut_, "\r\n");
}

void CCmdGauss9::printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                           int atInd1, int atInd2, int atInd3, int atInd4, double startAngle, double endAngle, double step) const
{
    C3DVector r;
    CAtom* atomPtr;
    int count;
    
    fprintf(stdOut_, "\t%%Chk = CheckPnt\r\n");
    fprintf(stdOut_, "\t# %s opt=ModRed nosymm\r\n\r\n\t%s\r\n\r\n\t%g %i\r\n", baseSetSpec.data(), heading.data(), charge, spinMul);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        if(atomPtr)
        {
            std::string ID = atomPtr->getID();
            if(state_->getCurrFrameIndex() < atomPtr->r_.size())
            {
                r = atomPtr->r_[state_->getCurrFrameIndex()];
                fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), r.x_, r.y_, r.z_);
            }
            else
            {
                fprintf(stdOut_, "\t%-5s%-13s%-13s%-13s\r\n", ID.data(), "na", "na", "na");
            }
        }
    }
    
    if(step != 0.0) count = (int)ceil((endAngle - startAngle) / step);
    else count = 1;

    fprintf(stdOut_, "\r\n\t%i %i %i %i %.4f S %i %.4f\r\n", atInd1+1, atInd2+1, atInd3+1, atInd4+1, startAngle, count, step);
}

void CCmdGauss9::printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                           int atInd1, int atInd2, int atInd3, double startAngle, double endAngle, double step) const
{
    C3DVector r;
    CAtom* atomPtr;
    int count;
    
    fprintf(stdOut_, "\t%%Chk = CheckPnt\r\n");
    fprintf(stdOut_, "\t# %s opt=ModRed nosymm\r\n\r\n\t%s\r\n\r\n\t%g %i\r\n", baseSetSpec.data(), heading.data(), charge, spinMul);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        if(atomPtr)
        {
            std::string ID = atomPtr->getID();
            if(state_->getCurrFrameIndex() < atomPtr->r_.size())
            {
                r = atomPtr->r_[state_->getCurrFrameIndex()];
                fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), r.x_, r.y_, r.z_);
            }
            else
            {
                fprintf(stdOut_, "\t%-5s%-13s%-13s%-13s\r\n", ID.data(), "na", "na", "na");
            }
        }
    }
    
    if(step != 0.0) count = (int)ceil((endAngle - startAngle) / step);
    else count = 1;
    
    fprintf(stdOut_, "\r\n\t%i %i %i %.4f S %i %.4f\r\n", atInd1+1, atInd2+1, atInd3+1, startAngle, count, step);
}

void CCmdGauss9::printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                           int atInd1, int atInd2, double startDist, double endDist, double step) const
{
    C3DVector r;
    CAtom* atomPtr;
    int count;
    
    fprintf(stdOut_, "\t%%Chk = CheckPnt\r\n");
    fprintf(stdOut_, "\t# %s opt=ModRed nosymm\r\n\r\n\t%s\r\n\r\n\t%g %i\r\n", baseSetSpec.data(), heading.data(), charge, spinMul);
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        if(atomPtr)
        {
            std::string ID = atomPtr->getID();
            if(state_->getCurrFrameIndex() < atomPtr->r_.size())
            {
                r = atomPtr->r_[state_->getCurrFrameIndex()];
                fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), r.x_, r.y_, r.z_);
            }
            else
            {
                fprintf(stdOut_, "\t%-5s%-13s%-13s%-13s\r\n", ID.data(), "na", "na", "na");
            }
        }
    }
    
    if(step != 0.0) count = (int)ceil((endDist - startDist) / step);
    else count = 1;
    
    if(count < 0) count = -count;
    
    fprintf(stdOut_, "\r\n\t%i %i %.4f S %i %.4f\r\n", atInd1+1, atInd2+1, startDist, count, step);
}

void CCmdGauss9::parseGenxyzfrominputCommand(std::string commandLine, int& arg)
{
    std::vector<std::shared_ptr<CAtom>> atomsList;
    std::string text;
    std::string line;
    bool lastLineInFile = false;
    bool lastLineEmpty = false;

    text = CASCIIUtility::getWord(commandLine, arg++);
    
    FILE* filePtr = fopen(text.data(), "r");
    if(filePtr)
    {
        int emptyLineCount = 1;
        int nonEmptyLineCount = 0;
        
        printf("\r\n");
        
        do
        {
            line = CFileUtility::readLine(filePtr, lastLineInFile);

            // When the number of empty lines read is equal to 3, then we
            // are at the coordinate section!
            if(emptyLineCount == 3)
            {
                if(!CASCIIUtility::isLineEmpty(line))
                {
                    // The first line of the coordinate section is the
                    // charge and spin multiplicity, so only start after
                    // the first line of this section!
                    if(nonEmptyLineCount > 0)
                    {
                        auto atomPtr = std::make_shared<CAtom>();
                        
                        text = CASCIIUtility::getWord(line, 0);
                        atomPtr->setID(text.data());
                        
                        text = CASCIIUtility::getWord(line, 1);
                        atomPtr->r_[0].x_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 2);
                        atomPtr->r_[0].y_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 3);
                        atomPtr->r_[0].z_ = atof(text.data());
                        
                        atomsList.emplace_back(atomPtr);
                    }
                    
                    nonEmptyLineCount++;
                }
            }
            
            // When 4 empty lines have been read, then we
            // know that a coordinate section has been read
            // and we should process it! Similarly, we know
            // that a final coordinate section has been read
            // if bLastLineInFile=true and aAtomsList.size()>0
            if(((emptyLineCount == 4) && !lastLineEmpty) || (lastLineInFile && (atomsList.size() > 0)))
            {
                fprintf(stdOut_, "\t%i\r\n\tXYZ file\r\n", (int)atomsList.size());
                for(int i=0; i<atomsList.size(); i++)
                {
                    text = atomsList[i]->getID();
                    fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", text.data(), atomsList[i]->r_[0].x_, atomsList[i]->r_[0].y_, atomsList[i]->r_[0].z_);
                }
                
                atomsList.clear();
                emptyLineCount = 0;
            }
            
            // Keep track of the number of encountered empty lines in the file,
            // and remember to count two or more consecutive empty lines as one!
            if(CASCIIUtility::isLineEmpty(line))
            {
                if(!lastLineEmpty) emptyLineCount++;
                lastLineEmpty = true;
                nonEmptyLineCount = 0;
            }
            else
            {
                lastLineEmpty = false;
            }
            
        } while(!lastLineInFile);
        
        fclose(filePtr);
    }
    else
    {
        printf("Error! Could not open file: %s!", text.data());
    }    
}

void CCmdGauss9::parseGenoptefromoutputCommand(std::string commandLine, int& arg)
{
    std::string text, line;
    bool lastLineInFile = false;
    bool foundOptCompl = false;
    double lastEnergy = 0.0;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    FILE* filePtr = fopen(text.data(), "r");
    if(filePtr)
    {
        printf("\r\n");
        
        do
        {
            line = CFileUtility::readLine(filePtr, lastLineInFile);
            
            if(line.find("SCF Done") != std::string::npos)
            {
                text = CASCIIUtility::getWord(line, 4);
                lastEnergy = atof(text.data());
                foundOptCompl = false;
            }
            if(line.find("Optimization completed.") != std::string::npos)
            {
                fprintf(stdOut_, "\tEOpt = %.10f a.u.\r\n", lastEnergy);
                foundOptCompl = true;
            }
            if(!foundOptCompl && (line.find("Normal termination") != std::string::npos))
            {
                fprintf(stdOut_, "\tEOpt = %.10f a.u.\r\n", lastEnergy);
            }
            
        } while(!lastLineInFile);
        
        fclose(filePtr);
    }
    else
    {
        printf("Error! Could not open file: %s!", text.data());
    }
}

void CCmdGauss9::parseGenoptxyzfromoutputCommand(std::string commandLine, int& arg)
{
    bool lastLineInFile = false;
    bool foundSymbZMatrix = false;
    bool foundOptCompl = false;
    bool end;
    std::string text;
    std::string line;
    C3DVector v;
    std::vector<std::string> IDs;
    std::vector<C3DVector> atoms;
    
    atoms.reserve(1000);
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    FILE* filePtr = fopen(text.data(), "r");
    if(filePtr)
    {
        printf("\r\n");
        
        do
        {
            line = CFileUtility::readLine(filePtr, lastLineInFile);
            
            if(foundSymbZMatrix && (line.find("Input orientation:") != std::string::npos))
            {
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);

                atoms.clear();
                end = false;
                do
                {
                    line = CFileUtility::readLine(filePtr, lastLineInFile);

                    if(line.find("------") == std::string::npos)
                    {
                        text = CASCIIUtility::getWord(line, 3);
                        v.x_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 4);
                        v.y_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 5);
                        v.z_ = atof(text.data());
                        
                        atoms.emplace_back(v);
                    }
                    else
                    {
                        end = true;
                    }
                    
                } while (!lastLineInFile && !end);

                foundOptCompl = false;
            }
            if(line.find("Symbolic Z-matrix:") != std::string::npos)
            {
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                
                IDs.clear();
                end = false;
                do
                {
                    line = CFileUtility::readLine(filePtr, lastLineInFile);

                    if(!CASCIIUtility::isLineEmpty(line))
                    {
                        text = CASCIIUtility::getWord(line, 0);
                        IDs.emplace_back(text);
                    }
                    else
                    {
                        end = true;
                    }

                } while(!lastLineInFile && !end);

                foundSymbZMatrix = true;
            }
            if(line.find("Optimization completed.") != std::string::npos)
            {
                printXYZ(IDs, atoms);
                foundOptCompl = true;
            }
            if(!foundOptCompl && (line.find("Normal termination") != std::string::npos))
            {
                printXYZ(IDs, atoms);
            }
            
        } while(!lastLineInFile);
        
        fclose(filePtr);
    }
    else
    {
        printf("Error! Could not open file: %s!", text.data());
    }
}

void CCmdGauss9::parseGenxyzfromoutputCommand(std::string commandLine, int& arg)
{
    bool lastLineInFile = false;
    bool foundSymbZMatrix = false;
    bool end;
    std::string text;
    std::string line;
    C3DVector v;
    std::vector<std::string> IDs;
    std::vector<C3DVector> atoms;
    
    atoms.reserve(1000);
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    FILE* filePtr = fopen(text.data(), "r");
    if(filePtr)
    {
        printf("\r\n");
        
        do
        {
            line = CFileUtility::readLine(filePtr, lastLineInFile);
            
            if(foundSymbZMatrix && (line.find("Input orientation:") != std::string::npos))
            {
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                
                atoms.clear();
                end = false;
                do
                {
                    line = CFileUtility::readLine(filePtr, lastLineInFile);
                    
                    if(line.find("------") == std::string::npos)
                    {
                        text = CASCIIUtility::getWord(line, 3);
                        v.x_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 4);
                        v.y_ = atof(text.data());
                        
                        text = CASCIIUtility::getWord(line, 5);
                        v.z_ = atof(text.data());
                        
                        atoms.emplace_back(v);
                    }
                    else
                    {
                        end = true;
                    }
                    
                } while (!lastLineInFile && !end);
                
                printXYZ(IDs, atoms);
            }
            if(line.find("Symbolic Z-matrix:") != std::string::npos)
            {
                line = CFileUtility::readLine(filePtr, lastLineInFile);
                
                IDs.clear();
                end = false;
                do
                {
                    line = CFileUtility::readLine(filePtr, lastLineInFile);
                    
                    if(!CASCIIUtility::isLineEmpty(line))
                    {
                        text = CASCIIUtility::getWord(line, 0);
                        IDs.emplace_back(text);
                    }
                    else
                    {
                        end = true;
                    }
                    
                } while(!lastLineInFile && !end);
                
                foundSymbZMatrix = true;
            }
            
        } while(!lastLineInFile);
        
        fclose(filePtr);
    }
    else
    {
        printf("Error! Could not open file: %s!", text.data());
    }
}

void CCmdGauss9::printXYZ(const std::vector<std::string> &IDs, const std::vector<C3DVector> &atoms) const
{
    fprintf(stdOut_, "\t%i\r\n\tXYZ file\r\n", (int)IDs.size());
    for(int i=0; i<IDs.size(); i++)
    {
        if(i < atoms.size())
        {
            fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", IDs[i].data(), atoms[i].x_, atoms[i].y_, atoms[i].z_);
        }
    }
}
