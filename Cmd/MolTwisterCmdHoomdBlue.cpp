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
#include <typeinfo>
#include "Utilities/LennardJonesUnits.h"
#include "Utilities/3DRect.h"
#include "MolTwisterCmdHoomdBlue.h"
#include "Tools/MolTwisterStateTools.h"
#include "MDFF/Bonds/MolTwisterMDFFBond_Harm.h"
#include "MDFF/Bonds/MolTwisterMDFFBond_LJC.h"
#include "MDFF/Bonds/MolTwisterMDFFBond_Morse.h"
#include "MDFF/Angles/MolTwisterMDFFAngle_Harm.h"
#include "MDFF/Dihedrals/MolTwisterMDFFDih_Harm.h"
#include "MDFF/Dihedrals/MolTwisterMDFFDih_Fourier4t.h"
#include "MDFF/Non-Bonded/MolTwisterMDFFNonBonded_Buck.h"
#include "MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ1208.h"
#include "MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJ.h"
#include "MDFF/Non-Bonded/MolTwisterMDFFNonBonded_LJBuck.h"

void CCmdHoomdblue::onAddKeywords()
{
    addKeyword("hoomdblue");
    addKeyword("genff");
    addKeyword("gendata");
    addKeyword("genrun");
    addKeyword("atmtolj");
    addKeyword("fstolj");
    addKeyword("Ktolj");
    addKeyword("bondacrosspbc");
    addKeyword("enable12and13");
}

std::string CCmdHoomdblue::getHelpString() const
{
    std::string text;

    text+= "\tUsage: hoomdblue genff\r\n";
    text+= "\t       hoomdblue gendata [bondacrosspbc]\r\n";
    text+= "\t       hoomdblue genrun <data file> <ff file> <temperature (K)> <time step (fs)> [enable12and13]\r\n";
    text+= "\t       hoomdblue atmtolj <value in atm units>\r\n";
    text+= "\t       hoomdblue fstolj <value in fs units>\r\n";
    text+= "\t       hoomdblue Ktolj <value in K units>\r\n";
    text+= "\r\n";
    text+= "\tThis command provides features that relate to version 1.3.1 of the HOOMD-blue MD simulator,\r\n";
    text+= "\tas well as other versions that support the same input / output structure. The 'atmtolj',\r\n";
    text+= "\t'fstolj' and 'Ktolj' sub-commands are used to convert values given in atm, fs and K,\r\n";
    text+= "\trespectiviely, to the Lennard-Jones units used by the sub-commands that generate inpu files\r\n";
    text+= "\tfor HOOMD-blue.\r\n";
    text+= "\r\n";
    text+= "\tThe 'genff' sub-command will create a HOOMD-blue Python text that sets up the specified force\r\n";
    text+= "\tfield. This should be piped to a Python file (i.e.,genff > <ff file>.py).\r\n";
    text+= "\r\n";
    text+= "\tThe 'gendata' sub-command will create an XML text that sets up all atomic positions, atomic IDs,\r\n";
    text+= "\tbond definitions, angle definitions, etc. This should be piped to a XML file (i.e.,gendata >\r\n";
    text+= "\t<data file>.xml). To make sure bonds are created across PBCs, use the 'bondacrosspbc' keyword.\r\n";
    text+= "\r\n";
    text+= "\tOnce 'genff' and 'gendata' has generated the <ff file>.py and <data file>.xml, the 'genrun'\r\n";
    text+= "\tsub-command can be applied to generate a Python run script for HOOMD-blue. This should be piped\r\n";
    text+= "\tto a Python file (i.e., genrun <data file> <ff file> <temp> <step> > run.py). To enable 1-2 and\r\n";
    text+= "\t1-3 interactions in HOOMD-blue, use the 'enable12and13' keyword.\r\n";
    text+= "\r\n";
    text+= "\tThe generated input files can now be used as a basis for creating the final input files.";

    return text;
}

void CCmdHoomdblue::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
     
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "genff")
    {
        parseGenffCommand(commandLine, arg);
    }
    
    else if(text == "gendata")
    {
        parseGendataCommand(commandLine, arg);
    }

    else if(text == "genrun")
    {
        parseGenrunCommand(commandLine, arg);
    }

    else if(text == "atmtolj")
    {
        parseAtmtoljCommand(commandLine, arg);
    }

    else if(text == "fstolj")
    {
        parseFstoljCommand(commandLine, arg);
    }

    else if(text == "Ktolj")
    {
        parseKtoljCommand(commandLine, arg);
    }
    
    else
    {
        printf("Syntax Error: Second argument should specify the kind of Hoomd-Blue script to generate!");
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1)  state_->view3D_->requestUpdate(true);
        else                                state_->view3D_->requestUpdate(false);
    }
}

void CCmdHoomdblue::parseGenrunCommand(std::string commandLine, int& arg)
{
    CLJUnits lju;
    std::string dataFileName, ffFileName, text;
    double timeStep, temperature;
    bool ignore12And13PairInt = true;

    
    dataFileName = CASCIIUtility::getWord(commandLine, arg++);
    ffFileName = CASCIIUtility::getWord(commandLine, arg++);

    text = CASCIIUtility::getWord(commandLine, arg++);
    temperature = atof(text.data());
    temperature/= lju.tempUnit();

    text = CASCIIUtility::getWord(commandLine, arg++);
    timeStep = atof(text.data());
    timeStep/= lju.timeUnit();
    
    
    // Collect a maximum of 10 optional specifications
    for(int i=0; i<10; i++)
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        if(text == "enable12and13") ignore12And13PairInt = false;
    }
    
    
    // Mandatory start of python script for Hoomd-Blue
    fprintf(stdOut_, "\tfrom hoomd_script import *\r\n");
    fprintf(stdOut_, "\tfrom hoomd_plugins import potentials\r\n");
    fprintf(stdOut_, "\tfrom math import sqrt\r\n");
    fprintf(stdOut_, "\timport random\r\n\r\n");
    
    
    // Read XML data file (generated by MolTwister gendata command)
    fprintf(stdOut_, "\t# Import data file (coordinates etc.)\r\n");
    fprintf(stdOut_, "\tcontext.initialize()\r\n");
    fprintf(stdOut_, "\tsystem = init.read_xml(filename=\"%s.xml\")\r\n\r\n", dataFileName.data());
    
    
    // Enable PPPM Ewald summation
    fprintf(stdOut_, "\t# Enable PPPM Ewald summation\r\n");
    fprintf(stdOut_, "\tcharged = group.charged()\r\n");
    fprintf(stdOut_, "\tpppm = charge.pppm(group=charged)\r\n");
    fprintf(stdOut_, "\tpppm.set_params(Nx=64, Ny=64, Nz=64, order=6, rcut=12.0)\r\n\r\n");

    
    // Import force-field description (generated by Hoomd-Blue genff command)
    fprintf(stdOut_, "\t# Import force-field file\r\n");
    fprintf(stdOut_, "\timport %s\r\n\r\n", ffFileName.data());
    
    
    // Construct a random velocity distribution at a temerature T
    fprintf(stdOut_, "\t# initialize the velocities to be a thermal distribution\r\n");
    fprintf(stdOut_, "\trandom.seed(1234);\r\n");
    fprintf(stdOut_, "\tT = %g\r\n", temperature);
    fprintf(stdOut_, "\tpx = py = pz = npart = 0.0;\r\n");
    fprintf(stdOut_, "\tfor p in system.particles:\r\n");
    fprintf(stdOut_, "\t\tmass = p.mass;\r\n");
    fprintf(stdOut_, "\t\tsigma = sqrt(T / mass)\r\n");
    fprintf(stdOut_, "\t\tvx = random.gauss(0, sigma)\r\n");
    fprintf(stdOut_, "\t\tvy = random.gauss(0, sigma)\r\n");
    fprintf(stdOut_, "\t\tvz = random.gauss(0, sigma)\r\n");
    fprintf(stdOut_, "\t\tp.velocity = (vx, vy, vz)\r\n\r\n");
    
    fprintf(stdOut_, "\t\tnpart = npart + 1.0;\r\n\r\n");
    
    fprintf(stdOut_, "\t\t# sum the total system momentum\r\n");
    fprintf(stdOut_, "\t\tpx += mass*vx;\r\n");
    fprintf(stdOut_, "\t\tpy += mass*vy;\r\n");
    fprintf(stdOut_, "\t\tpz += mass*vz;\r\n\r\n");
    
    fprintf(stdOut_, "\t# compute average momentum\r\n");
    fprintf(stdOut_, "\tpx /= npart;\r\n");
    fprintf(stdOut_, "\tpy /= npart;\r\n");
    fprintf(stdOut_, "\tpz /= npart;\r\n\r\n");
    
    fprintf(stdOut_, "\t# subtract that average momentum from each particle\r\n");
    fprintf(stdOut_, "\tfor p in system.particles:\r\n");
    fprintf(stdOut_, "\t\tmass = p.mass;\r\n");
    fprintf(stdOut_, "\t\tv = p.velocity;\r\n");
    fprintf(stdOut_, "\t\tp.velocity = (v[0] - px/mass, v[1] - py/mass, v[2] - pz/mass);\r\n\r\n");

    
    // Generate group of mobile atoms
    fprintf(stdOut_, "\t# Set up mobile region\r\n");
    std::vector<int> startTags, endTags;
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        bool bFoundMobiles = false;
        int iFirstTag = i;

        while((i<(int)state_->atoms_.size()) && state_->atoms_[i]->isMobile_)
        {
            bFoundMobiles = true;
            i++;
        }
        
        if(bFoundMobiles)
        {
            startTags.emplace_back(iFirstTag);
            endTags.emplace_back(i-1);
        }
    }
    
    for(int i=0; i<(int)startTags.size(); i++)
    {
        if(i == 0)
            fprintf(stdOut_, "\tranges=range(%i, %i)\r\n", startTags[i], endTags[i]+1);
        else
            fprintf(stdOut_, "\tranges=ranges+range(%i, %i)\r\n", startTags[i], endTags[i]+1);
    }
    fprintf(stdOut_, "\tunionMobile=group.tag_list(name=\"tag_list\", tags=ranges)\r\n\r\n");
    
    
    // Run simulation
    fprintf(stdOut_, "\t# Set up integration scheme and run simulation\r\n");
    if(ignore12And13PairInt) fprintf(stdOut_, "\tnlist.reset_exclusions(exclusions = ['1-2', '1-3', '1-4'])\r\n");
    fprintf(stdOut_, "\tintegrate.mode_standard(dt=%g)\r\n", timeStep);
    fprintf(stdOut_, "\tintegrate.nvt(group=unionMobile, T=%g, tau=%g)\r\n\r\n", temperature, 100.0*timeStep);
    
    fprintf(stdOut_, "\tdump.pdb(filename=\"initial.pdb\")\r\n");
    fprintf(stdOut_, "\tdump.xml(filename=\"final.xml\", restart=True, period=1000, all=True)\r\n");
    fprintf(stdOut_, "\tdump.dcd(filename=\"traj.dcd\", period=100)\r\n");
    fprintf(stdOut_, "\tlogdata=['time', 'temperature', 'pressure', 'potential_energy', 'kinetic_energy', 'volume']\r\n");
    fprintf(stdOut_, "\tlog=analyze.log(filename=\"output.log\", quantities=logdata, period=100, overwrite=True)\r\n\r\n");

    fprintf(stdOut_, "\tupdate.zero_momentum(period=100)\r\n\r\n");
    
    fprintf(stdOut_, "\trun(1e6)\r\n");
}

void CCmdHoomdblue::parseGenffCommand(std::string, int&)
{
    CMDFFNonBonded* nonBonded;
    CMDFFAngle* angle;
    CMDFFBond* bond;
    CMDFFDih* dih;
    std::vector<std::string> listOfAtomTypes;
    std::string stringAtom1, stringAtom2, stringAtom3, stringAtom4;
    bool havePrinted;
    bool hasLJ, hasLJ1208, hasBuck, hasMorse, hasHarmBond, hasLJCBond, hasHarmAngle, hasDihFour, hasDihHarm;
    
    
    state_->searchForAtomTypes(listOfAtomTypes);
    
    // Check enabled features
    hasLJ = hasLJ1208 = hasBuck = hasMorse = hasHarmBond = hasLJCBond = hasHarmAngle = hasDihFour = hasDihHarm = false;

    for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFNonBondedList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFNonBonded_LJ)) hasLJ = true;
        if(typeid(FFRef) == typeid(CMDFFNonBonded_LJBuck)) hasLJ = true;
    }

    for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFNonBondedList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFNonBonded_Buck)) hasBuck = true;
        if(typeid(FFRef) == typeid(CMDFFNonBonded_LJBuck)) hasBuck = true;
    }

    for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFNonBondedList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFNonBonded_LJ1208)) hasLJ1208 = true;
    }
    
    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFBondList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFBond_Morse)) hasMorse = true;
    }

    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFBondList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFBond_Harm)) hasHarmBond = true;
    }

    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFBondList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFBond_LJC)) hasLJCBond = true;
    }
    
    for(int i=0; i<state_->mdFFAngleList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFAngleList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFAngle_Harm)) hasHarmAngle = true;
    }

    for(int i=0; i<state_->mdFFDihList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFDihList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFDih_Fourier4t)) hasDihFour = true;
    }

    for(int i=0; i<state_->mdFFDihList_.size(); i++)
    {
        auto& FFRef = *state_->mdFFDihList_.get(i);
        if(typeid(FFRef) == typeid(CMDFFDih_Harm)) hasDihHarm = true;
    }
    

    // Import Hoomd-Blue libraries
    fprintf(stdOut_, "\tfrom hoomd_script import *\r\n");
    fprintf(stdOut_, "\tfrom hoomd_plugins import potentials\r\n\r\n");

    // Import python libraries
    fprintf(stdOut_, "\tfrom math import exp\r\n");
    fprintf(stdOut_, "\r\n");

    // Set up different definitions and pair styles
    fprintf(stdOut_, "\t# Set up different definitions and pair styles\r\n");
    if(hasLJ) fprintf(stdOut_, "\tlj = pair.lj(r_cut=12.0)\r\n");
    if(hasLJ1208) fprintf(stdOut_, "\tlj1208 = potentials.pair.lj1208(r_cut=12.0)\r\n");
    if(hasBuck) fprintf(stdOut_, "\tbuck = potentials.pair.buck(r_cut=12.0)\r\n");
    if(hasMorse) fprintf(stdOut_, "\tmorse = pair.morse(r_cut=12.0)\r\n");
    if(hasHarmBond) fprintf(stdOut_, "\tharmbond = bond.harmonic(name=\"harm_bond\")\r\n");
    if(hasLJCBond) fprintf(stdOut_, "\tljcbond = potentials.bond.lj_coul(name=\"ljc_bond\")\r\n");
    if(hasHarmAngle) fprintf(stdOut_, "\tharmangle = angle.harmonic()\r\n");
    if(hasDihFour) fprintf(stdOut_, "\tdihfour = dihedral.opls()\r\n");
    if(hasDihHarm) fprintf(stdOut_, "\tdihharm = dihedral.harmonic()\r\n");
    fprintf(stdOut_, "\r\n");
    
    // Print sorted pair coefficients list
    fprintf(stdOut_, "\t# Set all required pair coefficeints\r\n");
    havePrinted = false;
    for(int i=0; i<(int)listOfAtomTypes.size(); i++)
    {
        for(int j=i; j<(int)listOfAtomTypes.size(); j++)
        {
            if(hasMorse)
            {
                fprintf(stdOut_, "\tmorse.pair_coeff.set('%s', '%s', D0=0.0, alpha=0.0, r0=0.0)\r\n", listOfAtomTypes[i].data(), listOfAtomTypes[j].data());
                havePrinted = true;
            }
            if(hasLJ)
            {
                fprintf(stdOut_, "\tlj.pair_coeff.set('%s', '%s', epsilon=0.0, sigma=0.0)\r\n", listOfAtomTypes[i].data(), listOfAtomTypes[j].data());
                havePrinted = true;
            }
            if(hasLJ1208)
            {
                fprintf(stdOut_, "\tlj1208.pair_coeff.set('%s', '%s', epsilon=0.0, sigma=0.0)\r\n", listOfAtomTypes[i].data(), listOfAtomTypes[j].data());
                havePrinted = true;
            }
            if(hasBuck)
            {
                fprintf(stdOut_, "\tbuck.pair_coeff.set('%s', '%s', A=0.0, rho=1.0, C=0.0)\r\n", listOfAtomTypes[i].data(), listOfAtomTypes[j].data());
                havePrinted = true;
            }

            std::shared_ptr<std::vector<int>> indicesList = state_->mdFFNonBondedList_.indexFromNames(listOfAtomTypes[i], listOfAtomTypes[j]);
            for(int n=0; n<(int)indicesList->size(); n++)
            {
                nonBonded = state_->mdFFNonBondedList_.get((*indicesList)[n]);
                if(!nonBonded) continue;
                
                for(int k=0; k<nonBonded->getNumHoomdBlueDef(); k++)
                {
                    fprintf(stdOut_, "\t%s\r\n", nonBonded->getHoomdBlueDef(k).data());
                    havePrinted = true;
                }
            }
        }
    }

    // Print bond coefficients
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Declare bond coefficients to make sure all are defined\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        if(hasHarmBond)
        {
            fprintf(stdOut_, "\tharmbond.bond_coeff.set('bondtype%i', k=0.0, r0=0.0)\r\n", i+1);
            havePrinted = true;
        }
        if(hasLJCBond)
        {
            fprintf(stdOut_, "\tljcbond.bond_coeff.set('bondtype%i', epsilon=0.0, sigma=0.1, scale=0.0)\r\n", i+1);
            havePrinted = true;
        }
    }
    if(!havePrinted) fprintf(stdOut_, "\t# \tno declarations were needed...\r\n");

    int bondIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Set bond coefficients\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        bond = state_->mdFFBondList_.get(i);
        if(!bond) continue;
        
        bondIndex++;
        for(int k=0; k<bond->getNumHoomdBlueDef(); k++)
        {
            fprintf(stdOut_, "\t%s\r\n", bond->getHoomdBlueDef(k, bondIndex).data());
            havePrinted = true;
        }
    }

    // Print angle coefficients
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Declare angle coefficients to make sure all are defined\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFAngleList_.size(); i++)
    {
        if(hasHarmAngle)
        {
            fprintf(stdOut_, "\tharmangle.set_coeff('angletype%i', k=0.0, t0=0.0)\r\n", i+1);
            havePrinted = true;
        }
    }
    if(!havePrinted) fprintf(stdOut_, "\t# \tno declarations were needed...\r\n");

    int angleIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Set angle coefficients\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFAngleList_.size(); i++)
    {
        angle = state_->mdFFAngleList_.get(i);
        if(!angle) continue;
        
        angleIndex++;
        for(int k=0; k<angle->getNumHoomdBlueDef(); k++)
        {
            fprintf(stdOut_, "\t%s\r\n", angle->getHoomdBlueDef(k, angleIndex).data());
            havePrinted = true;
        }
    }
    
    // Print dihedral coefficients
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Declare dihedral coefficients to make sure all are defined\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFDihList_.size(); i++)
    {
        if(hasDihFour)
        {
            fprintf(stdOut_, "\tdihfour.set_coeff('dihedraltype%i', k1=0.0, k2=0.0, k3=0.0, k4=0.0)\r\n", i+1);
            havePrinted = true;
        }
        if(hasDihHarm)
        {
            fprintf(stdOut_, "\tdihharm.set_coeff('dihedraltype%i', k=0.0, d=1, n=0)\r\n", i+1);
            havePrinted = true;
        }
    }
    if(!havePrinted) fprintf(stdOut_, "\t# \tno declarations were needed...\r\n");

    int dihIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    fprintf(stdOut_, "\t# Set dihedral coefficients\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFDihList_.size(); i++)
    {
        dih = state_->mdFFDihList_.get(i);
        if(!dih) continue;
     
        dihIndex++;
        for(int k=0; k<dih->getNumHoomdBlueDef(); k++)
        {
            fprintf(stdOut_, "\t%s\r\n", dih->getHoomdBlueDef(k, dihIndex).data());
            havePrinted = true;
        }
    }
}

void CCmdHoomdblue::parseGendataCommand(std::string commandLine, int& arg)
{
    bool bondsAcrossPBC = false;
    std::string text;
    CAtom* atomPtr;
    C3DRect pbc = state_->view3D_->getPBC();
    CLJUnits lju;
    std::vector<int> bondAtoms1, bondAtoms2, bondMDTypeIndices;
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    std::vector<int> dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices;
    
    
    // Should we consider bonds across PBS?
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "bondacrosspbc") bondsAcrossPBC = true;
    else                                          arg--;
    
    
    // Find all bonds, angles and dihedrals
    CMolTwisterStateTools mtStateTools(state_, stdOut_);
    mtStateTools.getAllMDBondsInSystem(bondAtoms1, bondAtoms2, bondMDTypeIndices, bondsAcrossPBC);
    mtStateTools.getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices, bondsAcrossPBC);
    mtStateTools.getAllMDDihedralsInSystem(dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices, bondsAcrossPBC);
    
    
    // Print header
    printf("\r\n");
    fprintf(stdOut_, "\t<?xml version=\"1.0\" encoding=\"UTF-8\"?>\r\n");
    fprintf(stdOut_, "\t<hoomd_xml version=\"1.6\">\r\n");
    fprintf(stdOut_, "\t\t<configuration time_step=\"0\">\r\n");
    fprintf(stdOut_, "\t\t\t<box Lx=\"%.4f\" Ly=\"%.4f\" Lz=\"%.4f\"/>\r\n",
            pbc.getWidthX() / lju.distanceUnit(), pbc.getWidthY() / lju.distanceUnit(), pbc.getWidthZ() / lju.distanceUnit());
    
    
    // Print coordinates
    fprintf(stdOut_, "\t\t\t<position>\r\n");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        C3DVector pos;
        int iCurrFrame = state_->getCurrFrameIndex();
        atomPtr = state_->atoms_[i].get();
        if((iCurrFrame >= 0) && (iCurrFrame < (int)atomPtr->r_.size()))
            pos = atomPtr->r_[iCurrFrame];
        
        // [Pos] = AA, need input in m, so multimly with 10^{-10} before div. by dist units in m
        fprintf(stdOut_, "\t\t\t\t%.6f %.6f %.6f\r\n", pos.x_ / lju.distanceUnit(), pos.y_ / lju.distanceUnit(), pos.z_ / lju.distanceUnit());
    }
    fprintf(stdOut_, "\t\t\t</position>\r\n");
    
    
    // Print types
    fprintf(stdOut_, "\t\t\t<type>\r\n");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        std::string ID = atomPtr->getID();
        
        fprintf(stdOut_, "\t\t\t\t%s\r\n", ID.data());
    }
    fprintf(stdOut_, "\t\t\t</type>\r\n");
    
    
    // Print masses
    fprintf(stdOut_, "\t\t\t<mass>\r\n");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        
        // [m] = g/mol, need input in kg, so multimly with 10^{-3} and divide by NA before div. by mass units in kg/mol
        fprintf(stdOut_, "\t\t\t\t%.6f\r\n", atomPtr->m_ / lju.massUnit());
    }
    fprintf(stdOut_, "\t\t\t</mass>\r\n");

    
    // Print charges
    fprintf(stdOut_, "\t\t\t<charge>\r\n");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        atomPtr = state_->atoms_[i].get();
        
        // [Q] in units of |e|, need input in C, so multiply with |e| before div. by charge units in C
        fprintf(stdOut_, "\t\t\t\t%.6f\r\n", atomPtr->Q_ / lju.chargeUnit());
    }
    fprintf(stdOut_, "\t\t\t</charge>\r\n");
    
    
    // Print bonds
    if(bondAtoms1.size()) fprintf(stdOut_, "\t\t\t<bond>\r\n");
    for(int i=0; i<(int)bondAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t\t\t\tbondtype%i %i %i\r\n", bondMDTypeIndices[i]+1, bondAtoms1[i], bondAtoms2[i]);
    }
    if(bondAtoms1.size()) fprintf(stdOut_, "\t\t\t</bond>\r\n");

    
    // Print angles
    if(angleAtoms1.size()) fprintf(stdOut_, "\t\t\t<angle>\r\n");
    for(int i=0; i<(int)angleAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t\t\t\tangletype%i %i %i %i\r\n", angleMDTypeIndices[i]+1, angleAtoms1[i], angleAtoms2[i], angleAtoms3[i]);
    }
    if(angleAtoms1.size()) fprintf(stdOut_, "\t\t\t</angle>\r\n");

    
    // Print dihedrals
    if(dihAtoms1.size()) fprintf(stdOut_, "\t\t\t<dihedral>\r\n");
    for(int i=0; i<(int)dihAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t\t\t\tdihedraltype%i %i %i %i %i\r\n", dihMDTypeIndices[i]+1, dihAtoms1[i], dihAtoms2[i], dihAtoms3[i], dihAtoms4[i]);
    }
    if(dihAtoms1.size()) fprintf(stdOut_, "\t\t\t</dihedral>\r\n");

    
    // Print end
    fprintf(stdOut_, "\t\t</configuration>\r\n");
    fprintf(stdOut_, "\t</hoomd_xml>\r\n");
}

void CCmdHoomdblue::parseAtmtoljCommand(std::string commandLine, int& arg)
{
    CLJUnits lju;
    std::string text;
    double pressInAtm;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    pressInAtm = atof(text.data());
    
    printf("\r\n\t%g atm = %g in LJ like units\r\n", pressInAtm, pressInAtm / lju.pressUnit());
}

void CCmdHoomdblue::parseFstoljCommand(std::string commandLine, int& arg)
{
    CLJUnits lju;
    std::string text;
    double timeInFs;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    timeInFs = atof(text.data());
    
    printf("\r\n\t%g fs = %g in LJ like units\r\n", timeInFs, timeInFs / lju.timeUnit());
}

void CCmdHoomdblue::parseKtoljCommand(std::string commandLine, int& arg)
{
    CLJUnits lju;
    std::string text;
    double tempInK;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    tempInK = atof(text.data());
    
    printf("\r\n\t%g K = %g in LJ like units\r\n", tempInK, tempInK / lju.tempUnit());
}
