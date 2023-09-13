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
#include <vector>
#include "Utilities/3DVector.h"
#include "Utilities/3DRect.h"
#include "MolTwisterCmdPrint.h"
#include "Tools/MolTwisterStateTools.h"
#include "Tools/MolecularSystemTools.h"

void CCmdPrint::onAddKeywords()
{
    addKeyword("print");
    addKeyword("xyz");
    addKeyword("pdb");
    addKeyword("mtt");
    addKeyword("binary");
    addKeyword("bondsfromff");
    addKeyword("bondacrosspbc");
    addKeyword("currentframe");
    addKeyword("frame");
    addKeyword("version");
    addKeyword("mixff");
    addKeyword("resname");
    addKeyword("atomnames");
    addKeyword("atomicunits");
    addKeyword("tutorial");
}

std::string CCmdPrint::getHelpString() const
{
    std::string text;

    text+= "\tUsage: print <content>\r\n";
    text+= "\r\n";
    text+= "\tPrints, to screen, or file if piped ('>'), the specified contents, where <contents> can be any of the following\r\n";
    text+= "\t* xyz [currentframe] [atomicunits]: atomic IDs and positions in the XYZ-file format\r\n";
    text+= "\t* xyz [frame <frame index>] [atomicunits]: atomic IDs and positions in the XYZ-file format\r\n";
    text+= "\t* pdb [currentframe]: atomic IDs and positions in the PDB-file format (Protein Data Bank file format)\r\n";
    text+= "\t* pdb [frame <frame index>]: atomic IDs and positions in the PDB-file format (Protein Data Bank file format)\r\n";
    text+= "\t* mtt [bondsfromff] [bondacrosspbc] [currentframe] [binary]: atomic IDs and positions in the MolTwister file format\r\n";
    text+= "\t* mtt [bondsfromff] [bondacrosspbc] [frame <frame index>] [binary]: atomic IDs and positions in the MolTwister file format\r\n";
    text+= "\t* version: version of the software\r\n";
    text+= "\t* mixff resname <resname 1> resname <resname 2> <mixing-rule>: finds a mixed set of force field parameters between atoms in\r\n";
    text+= "\t                                                               <resname 1> and <resname 2> and outputs a list of mixed\r\n";
    text+= "\t                                                               parameters in the form of MolTwister exec() commands, that\r\n";
    text+= "\t                                                               constitutes a script that sets up the mixing force field\r\n";
    text+= "\t* mixff atomnames <list of atom IDs 1> atomnames <list of atom IDs 2> <mixing-rule>: finds a mixed set of force field parameters\r\n";
    text+= "\t                                                                                     between atoms in <list of atom IDs 1> and\r\n";
    text+= "\t                                                                                     <list of atom IDs 2> and outputs a list\r\n";
    text+= "\t                                                                                     of mixed parameters in the form of MolTwister\r\n";
    text+= "\t                                                                                     exec() commands, that constitutes a script\r\n";
    text+= "\t                                                                                     that sets up the mixing force field\r\n";
    text+= "\r\n";
    text+= "\tIf 'atomicunits' is not specified, then the positions are given in Aangstroms. If 'bondsfromff' is specified,\r\n";
    text+= "\tbonds are taken from the given force field bonds, otherwise bonds are calculated based on atomic separations\r\n";
    text+= "\tand will be identical to the bonds visible in the 3D view. If 'bondacrosspbc' is not specified, then the bonds\r\n";
    text+= "\twill not stretch across the periodic boundary conditions (PBC). The 'binary' keyword will invoke a binary\r\n";
    text+= "\tversion of the file format, if available. If neither 'currentframe' or 'frame' is specified, the default frame\r\n";
    text+= "\tto print is the frame at index 0. The <mixing-rule> parameter depends on the force fields to be mixed. In\r\n";
    text+= "\tsome cases it is not possible to choose the rule (the parameter will be ignored), but in others it is either\r\n";
    text+= "\t'aritmetic' or 'geometric'. Note that the <list of atom IDs X> consists of a comma separated list of atom IDs\r\n";
    text+= "\t(e.g., H, O, C7) without any space.";

    return text;
}

void CCmdPrint::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "xyz")
    {
        parseXYZCommand(commandLine, arg);
    }

    else if(text == "pdb")
    {
        parsePDBCommand(commandLine, arg);
    }

    else if(text == "mtt")
    {
        parseMTTCommand(commandLine, arg);
    }
    
    else if(text == "version")
    {
        parseVersionCommand(commandLine, arg);
    }

    else if(text == "mixff")
    {
        parseMixffCommand(commandLine, arg);
    }

    else
    {
        printf("Syntax Error: Second argument should specify what to print as!");
    }
}

void CCmdPrint::parseXYZCommand(std::string commandLine, int& arg)
{
    const double toAU = 1.8897261245651;
    C3DVector r;
    C3DRect pbc;
    std::string text;
    bool currentFrame = false;
    bool specificFrame = false;
    bool convertToAU = false;
    int frame = 0;

    if(state_->view3D_)
    {
        pbc = state_->view3D_->getPBC();
    }
    
    if(state_->atoms_.size() <= 0)
    {
        printf("\r\n\tNo atoms to print!\r\n");
        return;
    }

    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "currentframe")
    {
        currentFrame = true;
    }
    else if(text == "frame")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        frame = atoi(text.data());
        specificFrame = true;
        
    } else arg--;

    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "atomicunits")
    {
        convertToAU = true;
        
    } else arg--;
    
    printf("\r\n");
    for(int j=0; j<state_->atoms_[0]->r_.size(); j++)
    {
        if(currentFrame && (j != state_->currentFrame_)) continue;
        if(specificFrame && (j != frame)) continue;
        
        fprintf(stdOut_, "\t%i\r\n\tXYZ file, PBC = (%g, %g, %g)\r\n", (int)state_->atoms_.size(), pbc.getWidthX(), pbc.getWidthY(), pbc.getWidthZ());
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            std::string ID = state_->atoms_[i]->getID();
            if(j < state_->atoms_[i]->r_.size())
            {
                r = state_->atoms_[i]->r_[j];
                double x = r.x_;
                double y = r.y_;
                double z = r.z_;
                if(convertToAU) { x*= toAU; y*= toAU; z*= toAU;}
                fprintf(stdOut_, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), x, y, z);
            }
            else
            {
                fprintf(stdOut_, "\t%-5s%-13s%-13s%-13s\r\n", ID.data(), "na", "na", "na");
            }
        }
    }
    
    fprintf(stdOut_, "\r\n");
}

void CCmdPrint::parsePDBCommand(std::string commandLine, int& arg)
{
    C3DVector r;
    std::string text;
    bool currentFrame = false;
    bool specificFrame = false;
    char stringMolIndex[10] = "";
    char stringIndex[10] = "";
    int indexFrame = 0;
    
    if(state_->atoms_.size() <= 0)
    {
        printf("\r\n\tNo atoms to print!\r\n");
        return;
    }
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "currentframe")
    {
        currentFrame = true;
    }
    else if(text == "frame")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        indexFrame = atoi(text.data());
        specificFrame = true;
    }
    
    printf("\r\n");
    for(int j=0; j<state_->atoms_[0]->r_.size(); j++)
    {
        if(currentFrame && (j != state_->currentFrame_)) continue;
        if(specificFrame && (j != indexFrame)) continue;
        
        fprintf(stdOut_, "\tHEADER\r\n\tTITLE     Built with MolTwister\r\n");
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            int molIndex = state_->atoms_[i]->getMolIndex();
            if(molIndex < 10000)   sprintf(stringMolIndex, "%4i", molIndex);
            else                   sprintf(stringMolIndex, "****");
            
            if(i < 100000)  sprintf(stringIndex, "%5i", i);
            else            sprintf(stringIndex, "%5x", i);
            
            std::string ID = state_->atoms_[i]->getID();
            if(j < state_->atoms_[i]->r_.size())
            {
                r = state_->atoms_[i]->r_[j];
                fprintf(stdOut_, "\t%-6s%s %-4s %-3s A%s    %8.3f%8.3f%8.3f  0.00  0.00              \r\n", "ATOM",
                        stringIndex, ID.data(), state_->atoms_[i]->resname_.data(), stringMolIndex, r.x_, r.y_, r.z_);
            }
            else
            {
                r = state_->atoms_[i]->r_[j];
                fprintf(stdOut_, "\t%-6s%s %-4s %-3s A%s    %8s%8s%8s  0.00  0.00              \r\n", "ATOM",
                        stringIndex, ID.data(), state_->atoms_[i]->resname_.data(), stringMolIndex, "na", "na", "na");
            }
        }
    }
    
    fprintf(stdOut_, "\r\n");
}

void CCmdPrint::parseMTTCommand(std::string commandLine, int& arg)
{
    unsigned long ver = 1;
    C3DVector r;
    std::string text, resname;
    bool isCurrentFrame = false;
    bool isSpecificFrame = false;
    int indexFrame = 0;
    
    //////////////////////////////////////////////////////////////
    // *.mtt = MolTwister trajectory file
    //////////////////////////////////////////////////////////////
    // The file format is as follows. Line numbers
    // within the file are indicated by [n] followed
    // by a 'space', but are not part of the file.
    // -------------------------------------------------------
    // ASCII lines terminated by \r\n
    // -------------------------------------------------------
    // [1] 'mtt'
    // [2] <version number>
    // [3] 'binary' or 'ascii'
    // [4] '1' => molecular indices avaiable, '0' => not
    // [5] '1' => resnames avaiable, '0' => not
    // [6] '1' => bonds are avaiable, '0' => not
    // [7] '{'
    // [8] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    // [9] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    // [10]     .
    // [11]     .
    // [12] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    // [13] '}'
    // [14] '{'
    // [15]     .     next frame in trajectory...
    // [16]     .
    // [17] '}'
    // [18]     .
    // [19]     .
    // -------------------------------------------------------
    // Binary files (first three lines always ascii)
    //  note: string format is: <len(uchar)><...characters...>
    // -------------------------------------------------------
    // [1] 'mtt'
    // [2] <version number>
    // [3] 'binary' or 'ascii'
    // [4] 1 => molecular indices avaiable, 0 => not (1 byte)
    // [5] 1 => resnames avaiable, 0 => not (1 byte)
    // [6] 1 => bonds are avaiable, 0 => not (1 byte)
    // [7] -1 (long)
    // [8] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    //        (long)      (string)   (doubles)      (long)      (string)     (uchar)      (long)    ...     (long)
    // [9] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    //        (long)      (string)   (doubles)      (long)      (string)     (uchar)      (long)    ...     (long)
    // [10]     .
    // [11]     .
    // [12] <atom index> <atom type> <x> <y> <z> (<mol index>) (<resname>) (<#bonds N> <atom index1> ... <atom indexN>)
    //         (long)      (string)   (doubles)      (long)      (string)     (uchar)      (long)    ...     (long)
    // [13] -2 (long)
    // [14] -1 (long)
    // [15]     .     next frame in trajectory...
    // [16]     .
    // [17] -2 (long)
    // [18]     .
    // [19]     .
    //////////////////////////////////////////////////////////////
    // Note1! If resname is stated to be available, an empty
    //        resname is to be "---" in ASCII format. In binary
    //        format, the string length is simply zero
    // Note2! For ASCII files, empty <atom type> should always
    //        be written as "---"
    //////////////////////////////////////////////////////////////

    
    // Obtain command line parameters
    if(state_->atoms_.size() <= 0)
    {
        printf("\r\n\tNo atoms to print!\r\n");
        return;
    }

    bool bondInfoFromFF = false;
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "bondsfromff") bondInfoFromFF = true;
    else arg--;
    
    bool bondsAcrossPBC = false;
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "bondacrosspbc") bondsAcrossPBC = true;
    else arg--;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "currentframe")
    {
        isCurrentFrame = true;
    }
    else if(text == "frame")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        indexFrame = atoi(text.data());
        isSpecificFrame = true;
    } else arg--;

    bool isBinary = false;
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "binary") isBinary = true;
    else arg--;
    
    
    // Obtain bonds specified by MD force field (if required)
    std::vector<std::vector<int>> bondDestIndices;
    std::vector<int> molIndices;
    bondDestIndices.resize(state_->atoms_.size());
    if(bondInfoFromFF)
    {
        std::vector<int> mdBondsFromFF[2];
        std::vector<int> mdTypeIndex;
        CMolTwisterStateTools(state_, stdOut_).getAllMDBondsInSystem(mdBondsFromFF[0], mdBondsFromFF[1], mdTypeIndex, bondsAcrossPBC);
        for(int i=0; i<mdBondsFromFF[0].size(); i++)
        {
            if(mdBondsFromFF[0][i] < bondDestIndices.size())
                bondDestIndices[mdBondsFromFF[0][i]].emplace_back(mdBondsFromFF[1][i]);
            if(mdBondsFromFF[1][i] < bondDestIndices.size())
                bondDestIndices[mdBondsFromFF[1][i]].emplace_back(mdBondsFromFF[0][i]);
        }
        
        CMolecularSystemTools(state_, stdOut_).genMolIndices(bondDestIndices, molIndices);
    }
    
    
    // Print header
    printf("\r\n");
    if(!isBinary) fprintf(stdOut_, "\t");
    fprintf(stdOut_, "mtt\r\n");
    if(!isBinary) fprintf(stdOut_, "\t");
    fprintf(stdOut_, "%lu\r\n", ver);
    if(!isBinary) fprintf(stdOut_, "\t");
    fprintf(stdOut_, "%s\r\n", isBinary ? "binary" : "ascii");
    if(!isBinary)
    {
        fprintf(stdOut_, "\t%i\r\n", 1); // Molecular index available (0=false, 1=true)
        fprintf(stdOut_, "\t%i\r\n", 1); // Resname available (0=false, 1=true)
        fprintf(stdOut_, "\t%i\r\n", 1); // Bond info available (0=false, 1=true)
    }
    else
    {
        unsigned char byteWrite = 1;
        fwrite(&byteWrite, sizeof(byteWrite), 1, stdOut_); // Molecular index available (0=false, 1=true)
        fwrite(&byteWrite, sizeof(byteWrite), 1, stdOut_); // Resname available (0=false, 1=true)
        fwrite(&byteWrite, sizeof(byteWrite), 1, stdOut_); // Bond info available (0=false, 1=true)
    }
    
    
    // Print coordinate frames
    for(int j=0; j<state_->atoms_[0]->r_.size(); j++)
    {
        if(isCurrentFrame && (j != state_->currentFrame_)) continue;
        if(isSpecificFrame && (j != indexFrame)) continue;
     
        // Print start of coordinate frame
        if(!isBinary) fprintf(stdOut_, "\t{\r\n");
        else
        {
            long longWrite = -1;
            fwrite(&longWrite, sizeof(longWrite), 1, stdOut_); // Start of frame
        }

        // Print coordinate frame
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            std::string ID = state_->atoms_[i]->getID();
            resname = state_->atoms_[i]->resname_;
            long molIndex = -1;
            unsigned char numBonds;
            std::vector<long> bondsTo;
            if(!bondInfoFromFF)
            {
                molIndex = state_->atoms_[i]->getMolIndex();
                numBonds = (unsigned char)state_->atoms_[i]->getNumBonds();
                for(long l=0; l<numBonds; l++)
                {
                    CAtom* atomToBondTo = state_->atoms_[i]->getBondDest((int)l);
                    if(!atomToBondTo) continue;
                    
                    bondsTo.emplace_back(atomToBondTo->getAtomIndex());
                }
            }
            else
            {
                numBonds = (unsigned char)bondDestIndices[i].size();
                for(long l=0; l<numBonds; l++)
                    bondsTo.emplace_back(bondDestIndices[i][l]);
                if(i < molIndices.size()) molIndex = molIndices[i];
            }
            
            if(j < state_->atoms_[i]->r_.size())
            {
                r = state_->atoms_[i]->r_[j];
                
                if(!isBinary)
                {
                    text = resname;
                    CASCIIUtility::removeWhiteSpace(text);
                    if(text.size() == 0) resname = "---";

                    text = ID;
                    CASCIIUtility::removeWhiteSpace(text);
                    if(text.size() == 0) ID = "---";
                    
                    fprintf(stdOut_, "\t%li %s %g %g %g %li %s %i", (long)i, ID.data(), r.x_, r.y_, r.z_, molIndex, resname.data(), numBonds);
                    for(int l=0; l<bondsTo.size(); l++)
                        fprintf(stdOut_, " %li", bondsTo[l]);
                    fprintf(stdOut_, "\r\n");
                }
                else
                {
                    long l = (long)i;
                    fwrite(&l, sizeof(l), 1, stdOut_);                            // Atom index
                    unsigned char uc = (unsigned char)ID.size();
                    fwrite(&uc, sizeof(uc), 1, stdOut_);                          // String length
                    fwrite(ID.data(), sizeof(char), uc, stdOut_);                 // Atom type (O,C,Al, etc.)
                    fwrite(&r.x_, sizeof(r.x_), 1, stdOut_);                      // x-coordinate
                    fwrite(&r.y_, sizeof(r.y_), 1, stdOut_);                      // y-coordinate
                    fwrite(&r.z_, sizeof(r.z_), 1, stdOut_);                      // z-coordinate
                    fwrite(&molIndex, sizeof(molIndex), 1, stdOut_);              // Molecule index
                    uc = (unsigned char)resname.size();
                    fwrite(&uc, sizeof(uc), 1, stdOut_);                          // String length
                    fwrite(resname.data(), sizeof(char), uc, stdOut_);            // Residue name
                    fwrite(&numBonds, sizeof(numBonds), 1, stdOut_);              // Number of bonds
                    for(int l=0; l<bondsTo.size(); l++)
                        fwrite(&bondsTo[l], sizeof(bondsTo[l]), 1, stdOut_);      // Bond indices to connect to
                }
            }
        }

        // Print end of coordinate frame
        if(!isBinary) fprintf(stdOut_, "\t}\r\n");
        else
        {
            long longWrite = -2;
            fwrite(&longWrite, sizeof(longWrite), 1, stdOut_); // End of frame
        }
    }
    
    fflush(stdOut_);
    if(isBinary) skipTabRemoveReqFlag_ = true;
}

void CCmdPrint::parseVersionCommand(std::string, int&)
{
    fprintf(stdOut_, "\r\n\tVersion = %s\r\n", MOLTWISTER_VER);
}

void CCmdPrint::parseMixffCommand(std::string commandLine, int& arg)
{
    std::vector<std::string> constituents1, constituents2;
    std::string text, mixingRule;
    
    if(!getMixFFConstituents(commandLine, arg, constituents1)) return;
    if(!getMixFFConstituents(commandLine, arg, constituents2)) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text.size() == 0) { mixingRule = "none"; arg--; }
    else mixingRule = text;
    
    for(int i=0; i<constituents1.size(); i++)
    {
        std::shared_ptr<std::vector<int>> indicesI = state_->mdFFNonBondedList_.indexFromNames(constituents1[i], constituents1[i]);
        for(int m=0; m<indicesI->size(); m++)
        {
            CMDFFNonBonded* nonBondI = state_->mdFFNonBondedList_.get((*indicesI)[m]);
            if(!nonBondI) continue;
            
            for(int j=0; j<constituents2.size(); j++)
            {
                std::shared_ptr<std::vector<int>> indicesJ = state_->mdFFNonBondedList_.indexFromNames(constituents2[j], constituents2[j]);
                for(int n=0; n<indicesJ->size(); n++)
                {
                    CMDFFNonBonded* nonBondJ = state_->mdFFNonBondedList_.get((*indicesJ)[n]);
                    if(!nonBondJ) continue;
                    
                    std::string stringRet = nonBondI->molTwisterMixCmdArg(*nonBondJ, mixingRule);
                    
                    fprintf(stdOut_, "\texec(\"add mdnonbonded %s\");\r\n", stringRet.data());
                }
            }
        }
    }
}

bool CCmdPrint::getMixFFConstituents(std::string commandLine, int& arg, std::vector<std::string>& constituents) const
{
    std::vector<std::string> listOfAtomTypes;
    std::vector<std::string> listOfResnames;
    std::string text;

    constituents.clear();
    state_->searchForAtomTypes(listOfAtomTypes, &listOfResnames);
    
    printf("\r\n");
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "resname")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        for(int i=0; i<listOfResnames.size(); i++)
        {
            if(listOfResnames[i] == text)
            {
                constituents.emplace_back(listOfAtomTypes[i]);
            }
        }
    }
    else if(text == "atomnames")
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        CASCIIUtility::removeWhiteSpace(text);
        constituents = CASCIIUtility::getWords(text, ",");
    }
    else
    {
        printf("Error: expected 'resname' or 'atomnames'!");
        return false;
    }
    
    return true;
}
