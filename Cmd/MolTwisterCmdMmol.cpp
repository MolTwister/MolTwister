#include <iostream>
#include <vector>
#include <math.h>
#include "Utilities/3DVector.h"
#include "MolTwisterCmdMmol.h"

void CCmdMmol::onAddKeywords()
{
    addKeyword("mmol");
    addKeyword("toscript");
    addKeyword("noself");
    addKeyword("mixgeom");
}

std::string CCmdMmol::getHelpString() const
{
    std::string text;

    text+= "\tUsage: mmol toscript [noself] [mixgeom] <molfile 1> <molfile 2> ... <molfile N> > <script file to generate>\r\n";
    text+= "\r\n";
    text+= "\tThis command imports a series of *.mmol files, which are molecular definition\r\n";
    text+= "\tfiles defined within the MDynaMix MD software package, and converts them to a\r\n";
    text+= "\tseries of MolTwister exec(), script instructions that will set up the force\r\n";
    text+= "\tfields defined within the *.mmol files. The list of input molfiles (*.mmol)\r\n";
    text+= "\tis terminated by the pipe symbol, '>'.\r\n";
    text+= "\r\n";
    text+= "\tSelf interactions can be exculded in the MolTwister instructions by appyling the\r\n";
    text+= "\t'noself' keyword. By default, mixing of short range interaction parameters is done\r\n";
    text+= "\tby using arithmetic mixing rules. However, by employing the 'mixgeom' keyword,\r\n";
    text+= "\tgeometric mixing rules will be applied.";

    return text;
}

void CCmdMmol::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "toscript")
    {
        parseToscriptCommand(commandLine, arg);
    }
    
    else
    {
        printf("Syntax Error: Second argument should specify the mmol action!");
    }
}

bool CCmdMmol::doesDihExist(std::string atom1, std::string atom2, std::string atom3, std::string atom4, const std::vector<SDih>& dihedrals) const
{
    bool exists = false;
    
    for(int i=0; i<dihedrals.size(); i++)
    {
        if(atom1 == dihedrals[i].atom1_)
        {
            if(atom2 == dihedrals[i].atom2_)
            {
                if(atom3 == dihedrals[i].atom3_)
                {
                    if(atom4 == dihedrals[i].atom4_)
                    {
                        exists = true;
                        break;
                    }
                }
            }
        }
        
        if(atom1 == dihedrals[i].atom4_)
        {
            if(atom2 == dihedrals[i].atom3_)
            {
                if(atom3 == dihedrals[i].atom2_)
                {
                    if(atom4 == dihedrals[i].atom1_)
                    {
                        exists = true;
                        break;
                    }
                }
            }
        }
    }
    
    return exists;
}

bool CCmdMmol::doesAngleExist(std::string atom1, std::string atom2, std::string atom3, const std::vector<SAngle>& angles) const
{
    bool exists = false;
    
    for(int i=0; i<angles.size(); i++)
    {
        if(atom1 == angles[i].atom1_)
        {
            if(atom2 == angles[i].atom2_)
            {
                if(atom3 == angles[i].atom3_)
                {
                    exists = true;
                    break;
                }
            }
        }
        
        if(atom1 == angles[i].atom3_)
        {
            if(atom2 == angles[i].atom2_)
            {
                if(atom3 == angles[i].atom1_)
                {
                    exists = true;
                    break;
                }
            }
        }
    }
    
    return exists;
}

bool CCmdMmol::doesIntExist(std::string atom1, std::string atom2, const std::vector<SPairInt>& pairInt) const
{
    bool exists = false;
    
    for(int i=0; i<pairInt.size(); i++)
    {
        if(atom1 == pairInt[i].atom1_)
        {
            if(atom2 == pairInt[i].atom2_)
            {
                exists = true;
                break;
            }
        }

        if(atom1 == pairInt[i].atom2_)
        {
            if(atom2 == pairInt[i].atom1_)
            {
                exists = true;
                break;
            }
        }
    }
    
    return exists;
}

bool CCmdMmol::doesAtomExist(std::string atom, const std::vector<std::string>& atoms) const
{
    bool exists = false;
    
    for(int i=0; i<atoms.size(); i++)
    {
        if(atom == atoms[i])
        {
            exists = true;
            break;
        }
    }
    
    return exists;
}

void CCmdMmol::parseToscriptCommand(std::string commandLine, int& arg)
{
    std::vector<std::string> mmolFiles;
    std::shared_ptr<CMolDB> molDB = std::make_shared<CMolDB>();
    std::string text;
    bool includeSelfInteractions = true;
    bool mixGeom = false;
    int count = 0;
    int len;
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Get specifiers
    ////////////////////////////////////////////////////////////////////////////
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "noself") includeSelfInteractions = false;
    else arg--;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "mixgeom") mixGeom = true;
    else arg--;


    ////////////////////////////////////////////////////////////////////////////
    // Collect MMOL file names and open them
    ////////////////////////////////////////////////////////////////////////////
    bool pipe;
    do
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        len = (int)text.size();
        pipe = (text == ">") ? true : false;
        if((len > 0) && !pipe)
        {
            mmolFiles.emplace_back(text);
        }
        
        count++;
    
    } while((len != 0) && (count < 200) && (!pipe));
    
    for(int i=0; i<mmolFiles.size(); i++)
    {
        if(!molDB->addMoleculeType(mmolFiles[i].data()))
        {
            printf("Warning, could not open file: %s!", mmolFiles[i].data());
        }
    }
    

    ////////////////////////////////////////////////////////////////////////////
    // Generate self interactions
    ////////////////////////////////////////////////////////////////////////////
    if(includeSelfInteractions)
    {
        for(int i=0; i<molDB->getNumMoleculeTypes(); i++)
        {
            std::vector<std::string> atoms;
            std::vector<SPairInt> bonds;
            std::vector<SAngle> angles;
            std::vector<SDih> dihedrals;
            
            if(i != 0) fprintf(stdOut_, "\r\n\r\n\r\n");

            CMolDB::CMolecule* moleculePtr = molDB->getMoleculeType(i);
            if(!moleculePtr) continue;
            
            // Generate masses
            fprintf(stdOut_, "\r\n");
            atoms.clear();
            for(int a=0; a<moleculePtr->getNumAtoms(); a++)
            {
                CMolDB::CAtom* atomPtr = moleculePtr->getAtom(a);
                if(!atomPtr) continue;
                
                if(!doesAtomExist(atomPtr->name_, atoms))
                {
                    fprintf(stdOut_, "\texec(\"mod mass name %s to %.6f\");\r\n", atomPtr->name_.data(), atomPtr->mass_);
                    atoms.emplace_back(atomPtr->name_);
                }
            }
            
            // Generate charges
            fprintf(stdOut_, "\r\n");
            atoms.clear();
            for(int a=0; a<moleculePtr->getNumAtoms(); a++)
            {
                CMolDB::CAtom* atomPtr = moleculePtr->getAtom(a);
                if(!atomPtr) continue;
                
                if(!doesAtomExist(atomPtr->name_, atoms))
                {
                    fprintf(stdOut_, "\texec(\"mod charge name %s to %.6f\");\r\n", atomPtr->name_.data(), atomPtr->q_);
                    atoms.emplace_back(atomPtr->name_);
                }
            }
            
            // Special handling of MMOL fSPC is required (internal parameters of MDynaMix used for bonds and angles)
            if(moleculePtr->isFlexibleSPC())
            {
                
                fprintf(stdOut_, "\texec(\"add mdbond Ow Hw1 Morse %.6f %.6f %.6f\");\r\n", 426.702, 2.566, 1.0);
                fprintf(stdOut_, "\texec(\"add mdbond Ow Hw2 Morse %.6f %.6f %.6f\");\r\n", 426.702, 2.566, 1.0);
                fprintf(stdOut_, "\texec(\"add mdbond Hw1 Hw2 Harm %.6f %.6f\");\r\n", 687.41, 1.633);

                fprintf(stdOut_, "\texec(\"add mdangle Hw1 Hw2 Ow Class2 0.0 0.0 0.0 0.0 %.6f %.6f %.6f 0.0 0.0 0.0 0.0\");\r\n", -884.63, 1.633, 1.0);
                fprintf(stdOut_, "\texec(\"add mdangle Ow Hw1 Hw2 Class2 0.0 0.0 0.0 0.0 %.6f %.6f %.6f 0.0 0.0 0.0 0.0\");\r\n", -884.63, 1.633, 1.0);
                fprintf(stdOut_, "\texec(\"add mdangle Hw1 Ow Hw2 Class2 0.0 0.0 0.0 0.0 %.6f %.6f %.6f 0.0 0.0 0.0 0.0\");\r\n", 467.310, 1.0, 1.0);

                // Generate van der Waals interactions between atoms in molecule
                genVanDerWaals(moleculePtr, moleculePtr, mixGeom);
                continue;
            }
            
            // Generate bonds
            fprintf(stdOut_, "\r\n");
            bonds.clear();
            for(int b=0; b<moleculePtr->getNumBonds(); b++)
            {
                CMolDB::CBond* bondPtr = moleculePtr->getBond(b);
                if(!bondPtr) continue;
                
                int iConn1 = bondPtr->connection_[0] - 1;
                if((iConn1 < 0) || (iConn1 >= moleculePtr->getNumAtoms())) continue;

                int iConn2 = bondPtr->connection_[1] - 1;
                if((iConn2 < 0) || (iConn2 >= moleculePtr->getNumAtoms())) continue;
                
                CMolDB::CAtom* atom1Ptr = moleculePtr->getAtom(iConn1);
                if(!atom1Ptr) continue;

                CMolDB::CAtom* atom2Ptr = moleculePtr->getAtom(iConn2);
                if(!atom2Ptr) continue;
                
                if(!doesIntExist(atom1Ptr->name_, atom2Ptr->name_, bonds))
                {
                    if(bondPtr->type_ == 0) // Harmonic bond
                    {
                        fprintf(stdOut_, "\texec(\"add mdbond %s %s Harm %.6f %.6f\");\r\n", atom1Ptr->name_.data(), atom2Ptr->name_.data(), bondPtr->param2_, bondPtr->param1_);
                    }
                    else if(bondPtr->type_ == 1) // Morse bond
                    {
                        fprintf(stdOut_, "\texec(\"add mdbond %s %s Morse %.6f %.6f %.6f\");\r\n", atom1Ptr->name_.data(), atom2Ptr->name_.data(), bondPtr->param3_, bondPtr->param4_, bondPtr->param1_);
                    }
                    else
                    {
                        printf("Warning, MMOL type %i not supported (only Harmonic and Morse potentials are supported)!", bondPtr->type_);
                    }
                    
                    SPairInt paitInt(atom1Ptr->name_, atom2Ptr->name_);
                    bonds.emplace_back(paitInt);
                }
            }

            // Generate angles
            fprintf(stdOut_, "\r\n");
            angles.clear();
            for(int a=0; a<moleculePtr->getNumAngles(); a++)
            {
                CMolDB::CAngle* anglePtr = moleculePtr->getAngle(a);
                if(!anglePtr) continue;
                
                int iConn1 = anglePtr->connection_[0] - 1;
                if((iConn1 < 0) || (iConn1 >= moleculePtr->getNumAtoms())) continue;
                
                int iConn2 = anglePtr->connection_[1] - 1;
                if((iConn2 < 0) || (iConn2 >= moleculePtr->getNumAtoms())) continue;

                int iConn3 = anglePtr->connection_[2] - 1;
                if((iConn3 < 0) || (iConn3 >= moleculePtr->getNumAtoms())) continue;
                
                CMolDB::CAtom* atom1Ptr = moleculePtr->getAtom(iConn1);
                if(!atom1Ptr) continue;
                
                CMolDB::CAtom* atom2Ptr = moleculePtr->getAtom(iConn2);
                if(!atom2Ptr) continue;

                CMolDB::CAtom* atom3Ptr = moleculePtr->getAtom(iConn3);
                if(!atom3Ptr) continue;
                
                if(!doesAngleExist(atom1Ptr->name_, atom2Ptr->name_, atom3Ptr->name_, angles))
                {
                    fprintf(stdOut_, "\texec(\"add mdangle %s %s %s Harm %.6f %.6f\");\r\n", atom1Ptr->name_.data(), atom2Ptr->name_.data(), atom3Ptr->name_.data(), anglePtr->param2_, anglePtr->param1_);
                    
                    SAngle angle(atom1Ptr->name_, atom2Ptr->name_, atom3Ptr->name_);
                    angles.emplace_back(angle);
                }
            }

            // Generate dihedrals
            fprintf(stdOut_, "\r\n");
            dihedrals.clear();
            for(int d=0; d<moleculePtr->getNumDihedrals(); d++)
            {
                CMolDB::CDihedral* dihPtr = moleculePtr->getDihedral(d);
                if(!dihPtr) continue;
                
                int iConn1 = dihPtr->connection_[0] - 1;
                if((iConn1 < 0) || (iConn1 >= moleculePtr->getNumAtoms())) continue;
                
                int iConn2 = dihPtr->connection_[1] - 1;
                if((iConn2 < 0) || (iConn2 >= moleculePtr->getNumAtoms())) continue;
                
                int iConn3 = dihPtr->connection_[2] - 1;
                if((iConn3 < 0) || (iConn3 >= moleculePtr->getNumAtoms())) continue;

                int iConn4 = dihPtr->connection_[3] - 1;
                if((iConn4 < 0) || (iConn4 >= moleculePtr->getNumAtoms())) continue;
                
                CMolDB::CAtom* atom1Ptr = moleculePtr->getAtom(iConn1);
                if(!atom1Ptr) continue;
                
                CMolDB::CAtom* atom2Ptr = moleculePtr->getAtom(iConn2);
                if(!atom2Ptr) continue;
                
                CMolDB::CAtom* atom3Ptr = moleculePtr->getAtom(iConn3);
                if(!atom3Ptr) continue;

                CMolDB::CAtom* atom4Ptr = moleculePtr->getAtom(iConn4);
                if(!atom4Ptr) continue;
                
                if(!doesDihExist(atom1Ptr->name_, atom2Ptr->name_, atom3Ptr->name_, atom4Ptr->name_, dihedrals))
                {
                    fprintf(stdOut_, "\texec(\"add mddihedral %s %s %s %s Fourier4t %.6f %.6f %.6f %.6f\");\r\n",
                            atom1Ptr->name_.data(), atom2Ptr->name_.data(), atom3Ptr->name_.data(), atom4Ptr->name_.data(),
                            dihPtr->param1_, dihPtr->param2_, dihPtr->param3_, dihPtr->param4_);
                    
                    SDih dihedral(atom1Ptr->name_, atom2Ptr->name_, atom3Ptr->name_, atom4Ptr->name_);
                    dihedrals.emplace_back(dihedral);
                }
            }


            // Generate van der Waals interactions between atoms in molecule
            genVanDerWaals(moleculePtr, moleculePtr, mixGeom);
        }
    }
    

    ////////////////////////////////////////////////////////////////////////////
    // Generate molecular-molecular interactions
    ////////////////////////////////////////////////////////////////////////////
    for(int i=0; i<molDB->getNumMoleculeTypes(); i++)
    {
        for(int j=i+1; j<molDB->getNumMoleculeTypes(); j++)
        {
            CMolDB::CMolecule* molecule1Ptr = molDB->getMoleculeType(i);
            if(!molecule1Ptr) continue;
            CMolDB::CMolecule* molecule2Ptr = molDB->getMoleculeType(j);
            if(!molecule2Ptr) continue;

            genVanDerWaals(molecule1Ptr, molecule2Ptr, mixGeom);
        }
    }
}

void CCmdMmol::genVanDerWaals(const CMolDB::CMolecule* molecule1, const CMolDB::CMolecule* molecule2, bool mixGeom) const
{
    // Generate van der Waals (within different molecules, do not duplicate interactions with same atom names)
    // Buckinham potentials are prioritized over Lennard-Jones potentials, we chaeck if Buckingham
    // potentials are available by the non-zero entry of A and C in *.mmol files
    std::vector<SPairInt> pairInt;
    fprintf(stdOut_, "\r\n");
    for(int a1=0; a1<molecule1->getNumAtoms(); a1++)
    {
        const CMolDB::CAtom* atom1Ptr = molecule1->getAtom(a1);
        if(!atom1Ptr) continue;
        
        for(int a2=0; a2<molecule2->getNumAtoms(); a2++)
        {
            const CMolDB::CAtom* atom2Ptr = molecule2->getAtom(a2);
            if(!atom2Ptr) continue;
            
            if(!doesIntExist(atom1Ptr->name_, atom2Ptr->name_, pairInt))
            {
                bool isBuck = false;
                
                SPairInt pairIntInner(atom1Ptr->name_, atom2Ptr->name_);
                pairInt.emplace_back(pairIntInner);
                
                if((atom1Ptr->buckA_ != 0.0) || (atom1Ptr->buckC_ != 0.0))
                {
                    if((atom2Ptr->buckA_ != 0.0) || (atom2Ptr->buckC_ != 0.0))
                    {
                        isBuck = true;
                    }
                }
                
                if(isBuck)
                {
                    double A = sqrt(atom1Ptr->buckA_ * atom2Ptr->buckA_);
                    double rho = (atom1Ptr->buckRho_ + atom2Ptr->buckRho_) / 2.0;
                    double C = sqrt(atom1Ptr->buckC_ * atom2Ptr->buckC_);
                    fprintf(stdOut_, "\texec(\"add mdnonbonded %s %s Buck %.6f %.6f %.6f\");\r\n", atom1Ptr->name_.data(), atom2Ptr->name_.data(), A, rho, C);
                }
                else
                {
                    double epsilon;
                    double sigma;
                    if(mixGeom)
                    {
                        epsilon = sqrt(atom1Ptr->ljEpsilon_ * atom2Ptr->ljEpsilon_);
                        sigma = sqrt(atom1Ptr->ljSigma_ * atom2Ptr->ljSigma_);
                    }
                    else
                    {
                        epsilon = sqrt(atom1Ptr->ljEpsilon_ * atom2Ptr->ljEpsilon_);
                        sigma = (atom1Ptr->ljSigma_ + atom2Ptr->ljSigma_) / 2.0;
                    }
                    fprintf(stdOut_, "\texec(\"add mdnonbonded %s %s LJ %.6f %.6f\");\r\n", atom1Ptr->name_.data(), atom2Ptr->name_.data(), epsilon, sigma);
                }
            }
        }
    }
}
