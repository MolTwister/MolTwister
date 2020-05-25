#include <iostream>
#include <math.h>
#include <vector>
#include <string>
#include "Utilities/3DRect.h"
#include "MolTwisterCmdLammps.h"
#include "Tools/MolTwisterStateTools.h"

void CCmdLammps::onAddKeywords()
{
    addKeyword("lammps");
    addKeyword("genff");
    addKeyword("gendata");
    addKeyword("bondacrosspbc");
}

std::string CCmdLammps::getHelpString() const
{
    std::string text;

    text+= "\tUsage: lammps genff\r\n";
    text+= "\t       lammps gendata [bondacrosspbc]\r\n";
    text+= "\r\n";
    text+= "\tThis command provides features that relate to version 5 Sep 2014 version of the LAMMPS MD\r\n";
    text+= "\tsimulator, as well as other versions that support the same input / output structure.\r\n";
    text+= "\r\n";
    text+= "\tThe 'genff' sub-command will create a LAMMPS text that sets up the specified force field.\r\n";
    text+= "\tThis should be piped to a file (i.e.,genff > <ff file>).\r\n";
    text+= "\r\n";
    text+= "\tThe 'gendata' sub-command will create a LAMMPS text that sets up all atomic positions, atomic IDs,\r\n";
    text+= "\tbond definitions, angle definitions, etc. This should be piped to a file (i.e.,gendata >\r\n";
    text+= "\t<data file>). To make sure bonds are created across PBCs, use the 'bondacrosspbc' keyword.\r\n";
    text+= "\r\n";
    text+= "\tThe generated input files can now be used as a basis for creating the final input files. The\r\n";
    text+= "\t<ff file> and <data file> should be used as input to the main input file for LAMMPS.";

    return text;
}

void CCmdLammps::execute(std::string commandLine)
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

    else
    {
        printf("Syntax Error: Second argument should specify the kind of Lammps script to generate!");
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1)  state_->view3D_->requestUpdate(true);
        else                                state_->view3D_->requestUpdate(false);
    }
}

void CCmdLammps::parseGenffCommand(std::string, int&)
{
    CMDFFNonBonded* nonBonded;
    CMDFFAngle* angle;
    CMDFFBond* bond;
    CMDFFDih* dih;
    std::vector<std::string> listOfAtomTypes;
    std::string stringAtom1, stringAtom2, stringAtom3, stringAtom4;
    bool havePrinted;
    
    state_->searchForAtomTypes(listOfAtomTypes);

    // Print sorted pair coefficients list
    fprintf(stdOut_, "\r\n");
    havePrinted = false;
    for(int i=0; i<(int)listOfAtomTypes.size(); i++)
    {
        for(int j=i; j<(int)listOfAtomTypes.size(); j++)
        {
            std::shared_ptr<std::vector<int>> indicesList = state_->mdFFNonBondedList_.indexFromNames(listOfAtomTypes[i], listOfAtomTypes[j]);
            if(indicesList->size())
            {
                for(int n=0; n<(int)indicesList->size(); n++)
                {
                    nonBonded = state_->mdFFNonBondedList_.get((*indicesList)[n]);
                    if(!nonBonded) continue;
                    
                    for(int k=0; k<nonBonded->getNumLammpsDef(); k++)
                    {
                        fprintf(stdOut_, "\tpair_coeff %i %i %s\r\n", i+1, j+1, nonBonded->getLammpsDef(k, true).data());
                        havePrinted = true;
                    }
                }
            }
            else
            {
                fprintf(stdOut_, "\tpair_coeff %i %i lj/cut/coul/long 0.0 0.0     # %s %s\r\n", i+1, j+1, listOfAtomTypes[i].data(), listOfAtomTypes[j].data());
                havePrinted = true;
            }
        }
    }

    // Print bond coefficients
    int bondIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFBondList_.size(); i++)
    {
        bond = state_->mdFFBondList_.get(i);
        if(!bond) continue;

        bondIndex++;
        for(int k=0; k<bond->getNumLammpsDef(); k++)
        {
            fprintf(stdOut_, "\tbond_coeff %i %s\r\n", bondIndex, bond->getLammpsDef(k, true).data());
            havePrinted = true;
        }
    }

    // Print angle coefficients
    int angleIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFAngleList_.size(); i++)
    {
        angle = state_->mdFFAngleList_.get(i);
        if(!angle) continue;
        
        angleIndex++;
        for(int k=0; k<angle->getNumLammpsDef(); k++)
        {
            fprintf(stdOut_, "\tangle_coeff %i %s\r\n", angleIndex, angle->getLammpsDef(k, true).data());
            havePrinted = true;
        }
    }

    // Print dihedral coefficients
    int dihIndex = 0;
    if(havePrinted) fprintf(stdOut_, "\r\n\r\n");
    havePrinted = false;
    for(int i=0; i<state_->mdFFDihList_.size(); i++)
    {
        dih = state_->mdFFDihList_.get(i);
        if(!dih) continue;
        
        dihIndex++;
        for(int k=0; k<dih->getNumLammpsDef(); k++)
        {
            fprintf(stdOut_, "\tdihedral_coeff %i %s\r\n", dihIndex, dih->getLammpsDef(k, true).data());
            havePrinted = true;
        }
    }
}

void CCmdLammps::parseGendataCommand(std::string commandLine, int& arg)
{
    bool bondsAcrossPBC = false;
    std::string text;
    CAtom* atomPtr;
    C3DRect pbc = state_->view3D_->getPBC();
    std::vector<int> bondAtoms1, bondAtoms2, bondMDTypeIndices;
    std::vector<int> angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices;
    std::vector<int> dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices;
    std::vector<std::string> listOfAtomTypes;
    
    // Should we consider bonds across PBS?
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "bondacrosspbc") bondsAcrossPBC = true;
    else                        arg--;
    
    // Find all bonds, angles and dihedrals
    CMolTwisterStateTools mtStateTools(state_, stdOut_);
    mtStateTools.getAllMDBondsInSystem(bondAtoms1, bondAtoms2, bondMDTypeIndices, bondsAcrossPBC);
    mtStateTools.getAllMDAnglesInSystem(angleAtoms1, angleAtoms2, angleAtoms3, angleMDTypeIndices, bondsAcrossPBC);
    mtStateTools.getAllMDDihedralsInSystem(dihAtoms1, dihAtoms2, dihAtoms3, dihAtoms4, dihMDTypeIndices, bondsAcrossPBC);

    // Find all atom types
    state_->searchForAtomTypes(listOfAtomTypes);
    
    // Print header
    fprintf(stdOut_, "\r\n\tLAMMPS data file. CGCMM style. atom_style full generated by MolTwister\r\n");
    
    // Print num atoms, bonds, angles, dihedrals and impropers
    fprintf(stdOut_, "\t %i atoms\r\n", (int)state_->atoms_.size());
    fprintf(stdOut_, "\t %i bonds\r\n", (int)bondAtoms1.size());
    fprintf(stdOut_, "\t %i angles\r\n", (int)angleAtoms1.size());
    fprintf(stdOut_, "\t %i dihedrals\r\n", (int)dihAtoms1.size());
    fprintf(stdOut_, "\t %i impropers\r\n", 0);

    // Print num atom, bond, angle, dihedral and improper types
    fprintf(stdOut_, "\t %i atom types\r\n", (int)listOfAtomTypes.size());
    fprintf(stdOut_, "\t %i bond types\r\n", state_->mdFFBondList_.size());
    fprintf(stdOut_, "\t %i angle types\r\n", state_->mdFFAngleList_.size());
    fprintf(stdOut_, "\t %i dihedral types\r\n", state_->mdFFDihList_.size());
    fprintf(stdOut_, "\t %i improper types\r\n", 0);
    fprintf(stdOut_, "\t %.4f %.4f xlo xhi\r\n", pbc.rLow_.x_, pbc.rHigh_.x_);
    fprintf(stdOut_, "\t %.4f %.4f ylo yhi\r\n", pbc.rLow_.y_, pbc.rHigh_.y_);
    fprintf(stdOut_, "\t %.4f %.4f zlo zhi\r\n\r\n", pbc.rLow_.z_, pbc.rHigh_.z_);
    
    // Print masses
    fprintf(stdOut_, "\t Masses\r\n\r\n");
    for(int i=0; i<(int)listOfAtomTypes.size(); i++)
    {
        atomPtr = state_->getFirstOccurenceOf(listOfAtomTypes[i]);
        if(atomPtr)
        {
            fprintf(stdOut_, "\t%i %.6f     # %s\r\n", i+1, atomPtr->m_, listOfAtomTypes[i].data());
        }
    }

    // Print atoms
    fprintf(stdOut_, "\r\n\r\n\t Atoms\r\n\r\n");
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        C3DVector pos;
        int currFrame = state_->getCurrFrameIndex();
        atomPtr = state_->atoms_[i].get();
        std::string ID = atomPtr->getID();
        if((currFrame >= 0) && (currFrame < (int)atomPtr->r_.size()))
           pos = atomPtr->r_[currFrame];

        fprintf(stdOut_, "\t%i %i %i %.6f %.6f %.6f %.6f\t# %s\r\n",
                i+1, atomPtr->getMolIndex()+1, CMolTwisterState::atomTypeToTypeIndex(listOfAtomTypes, ID)+1, atomPtr->Q_, pos.x_, pos.y_, pos.z_, ID.data());
    }

    // Print bonds
    if(bondAtoms1.size()) fprintf(stdOut_, "\r\n\r\n\t Bonds\r\n\r\n");
    for(int i=0; i<(int)bondAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t%i %i %i %i\r\n", i+1, bondMDTypeIndices[i]+1, bondAtoms1[i]+1, bondAtoms2[i]+1);
    }

    // Print angles
    if(angleAtoms1.size()) fprintf(stdOut_, "\r\n\r\n\t Angles\r\n\r\n");
    for(int i=0; i<(int)angleAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t%i %i %i %i %i\r\n", i+1, angleMDTypeIndices[i]+1, angleAtoms1[i]+1, angleAtoms2[i]+1, angleAtoms3[i]+1);
    }

    // Print Dihedrals
    if(dihAtoms1.size()) fprintf(stdOut_, "\r\n\r\n\t Dihedrals\r\n\r\n");
    for(int i=0; i<(int)dihAtoms1.size(); i++)
    {
        fprintf(stdOut_, "\t%i %i %i %i %i %i\r\n", i+1, dihMDTypeIndices[i]+1, dihAtoms1[i]+1, dihAtoms2[i]+1, dihAtoms3[i]+1, dihAtoms4[i]+1);
    }
}
