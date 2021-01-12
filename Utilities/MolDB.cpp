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

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "Utilities/ASCIIUtility.h"
#include "Utilities/FileUtility.h"
#include "MolDB.h"

#define WHITESPACES  " \t,:"


////////////////////////////////////////////////////////////////////////
// CMolDB::CAtom
////////////////////////////////////////////////////////////////////////

CMolDB::CAtom::CAtom()
{
    x_ = 0.0;
    y_ = 0.0;
    z_ = 0.0;
    mass_ = 0.0;
    q_ = 0.0;
    ljSigma_ = 0.0;
    ljEpsilon_ = 0.0;
    buckA_ = 0.0;
    buckRho_ = 0.0;
    buckC_ = 0.0;
    
    atomType_ = -1;
    atomTypeGlobal_ = -1;
    firstAtomOfThisType_ = true;
}

bool CMolDB::CAtom::parse(std::string line)
{
    std::string word;
    
    // Retrieve atom name (e.g. C14)
    name_ = CASCIIUtility::getWord(line, 0, WHITESPACES);
    
    // Retrieve X, Y and Z initial local positions
    word = CASCIIUtility::getWord(line, 1, WHITESPACES);
    x_ = atof(word.data());
    word = CASCIIUtility::getWord(line, 2, WHITESPACES);
    y_ = atof(word.data());
    word = CASCIIUtility::getWord(line, 3, WHITESPACES);
    z_ = atof(word.data());
    
    // Retrieve atom mass
    word = CASCIIUtility::getWord(line, 4, WHITESPACES);
    mass_ = atof(word.data());

    // Retrieve atom charge
    word = CASCIIUtility::getWord(line, 5, WHITESPACES);
    q_ = atof(word.data());

    // Retrieve Lennard-Jones parameters
    word = CASCIIUtility::getWord(line, 6, WHITESPACES);
    ljSigma_ = atof(word.data());
    word = CASCIIUtility::getWord(line, 7, WHITESPACES);
    ljEpsilon_ = atof(word.data());

    // Retrieve Buckinham parameters, if available
    word = CASCIIUtility::getWord(line, 8, WHITESPACES);
    if(word[0] != '\0') buckA_ = atof(word.data());
    word = CASCIIUtility::getWord(line, 9, WHITESPACES);
    if(word[0] != '\0') buckRho_ = atof(word.data());
    word = CASCIIUtility::getWord(line, 10, WHITESPACES);
    if(word[0] != '\0') buckC_ = atof(word.data());
    
    return true;
}

bool CMolDB::CAtom::operator==(const CAtom& src) const
{
    if(fabs(mass_ - src.mass_) > 1E-9) return false;
    if(fabs(q_ - src.q_) > 1E-9) return false;
    if(fabs(ljSigma_ - src.ljSigma_) > 1E-9) return false;
    if(fabs(ljEpsilon_ - src.ljEpsilon_) > 1E-9) return false;
    if(fabs(buckA_ - src.buckA_) > 1E-9) return false;
    if(fabs(buckRho_ - src.buckRho_) > 1E-9) return false;
    if(fabs(buckC_ - src.buckC_) > 1E-9) return false;
    
    return true;
}

void CMolDB::CAtom::print() const
{
    printf("%s     %g %g %g     %g %g     %g %g     %g %g %g\r\n", 
           name_.data(), x_, y_, z_, mass_, q_,
           ljSigma_, ljEpsilon_, buckA_, buckRho_, buckC_);
}


////////////////////////////////////////////////////////////////////////
// CMolDB::CBond
////////////////////////////////////////////////////////////////////////

CMolDB::CBond::CBond()
{
    connection_[0] = -1;
    connection_[1] = -1;
    
    type_ = -1;
    param1_ = 0.0;
    param2_ = 0.0;
    param3_ = 0.0;
    param4_ = 0.0;
}

bool CMolDB::CBond::parse(std::string line)
{
    std::string word;
    
    // Retrieve bond type
    word = CASCIIUtility::getWord(line, 0, WHITESPACES);
    type_ = (char)atoi(word.data());
    
    // Retrieve connections
    word = CASCIIUtility::getWord(line, 1, WHITESPACES);
    connection_[0] = (char)atoi(word.data());
    
    word = CASCIIUtility::getWord(line, 2, WHITESPACES);
    connection_[1] = (char)atoi(word.data());
    
    // Retrieve parameters
    word = CASCIIUtility::getWord(line, 3, WHITESPACES);
    param1_ = atof(word.data());
    
    word = CASCIIUtility::getWord(line, 4, WHITESPACES);
    param2_ = atof(word.data());
    
    if(type_ == 1)
    {
        word = CASCIIUtility::getWord(line, 5, WHITESPACES);
        param3_ = atof(word.data());
        
        word = CASCIIUtility::getWord(line, 6, WHITESPACES);
        param4_ = atof(word.data());
    }
    
    return true;
}


////////////////////////////////////////////////////////////////////////
// CMolDB::CAngle
////////////////////////////////////////////////////////////////////////

CMolDB::CAngle::CAngle()
{
    connection_[0] = -1;
    connection_[1] = -1;
    connection_[2] = -1;
    
    param1_ = 0.0;
    param2_ = 0.0;
}

bool CMolDB::CAngle::parse(std::string line)
{
    std::string word;
        
    // Retrieve connections
    word = CASCIIUtility::getWord(line, 0, WHITESPACES);
    connection_[0] = (char)atoi(word.data());
    
    word = CASCIIUtility::getWord(line, 1, WHITESPACES);
    connection_[1] = (char)atoi(word.data());

    word = CASCIIUtility::getWord(line, 2, WHITESPACES);
    connection_[2] = (char)atoi(word.data());
    
    // Retrieve parameters
    word = CASCIIUtility::getWord(line, 3, WHITESPACES);
    param1_ = atof(word.data());
    
    word = CASCIIUtility::getWord(line, 4, WHITESPACES);
    param2_ = atof(word.data());
    
    return true;
}


////////////////////////////////////////////////////////////////////////
// CMolDB::CDihedral
////////////////////////////////////////////////////////////////////////

CMolDB::CDihedral::CDihedral()
{
    connection_[0] = -1;
    connection_[1] = -1;
    connection_[2] = -1;
    connection_[3] = -1;
    
    param1_ = 0.0;
    param2_ = 0.0;
    param3_ = 0.0;
    param4_ = 0.0;
    param5_ = 0.0;
    
    torsType_ = -1;
}

bool CMolDB::CDihedral::parse(std::string line, int torsType)
{
    std::string word;
        
    torsType_ = torsType;
    
    // Retrieve connections
    word = CASCIIUtility::getWord(line, 0, WHITESPACES);
    connection_[0] = (char)atoi(word.data());
    
    word = CASCIIUtility::getWord(line, 1, WHITESPACES);
    connection_[1] = (char)atoi(word.data());
    
    word = CASCIIUtility::getWord(line, 2, WHITESPACES);
    connection_[2] = (char)atoi(word.data());

    word = CASCIIUtility::getWord(line, 3, WHITESPACES);
    connection_[3] = (char)atoi(word.data());
    
    // Retrieve parameters
    word = CASCIIUtility::getWord(line, 4, WHITESPACES);
    param1_ = atof(word.data());
    
    word = CASCIIUtility::getWord(line, 5, WHITESPACES);
    param2_ = atof(word.data());

    word = CASCIIUtility::getWord(line, 6, WHITESPACES);
    param3_ = atof(word.data());
    
    if(torsType == 5)
    {
        word = CASCIIUtility::getWord(line, 7, WHITESPACES);
        param4_ = atof(word.data());

        word = CASCIIUtility::getWord(line, 8, WHITESPACES);
        param5_ = atof(word.data());
    }
    
    return true;
}


////////////////////////////////////////////////////////////////////////
// CMolDB::CMolecule
////////////////////////////////////////////////////////////////////////

bool CMolDB::CMolecule::load(FILE* fileHandle, int offsetAtomTypeGlobal)
{
    std::string line;
    bool lastLineInFile;
    int atomType, index;
    
    // Retrieve atoms
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);

    } while((line[0] == '#') && !lastLineInFile);

    CASCIIUtility::removeWhiteSpace(line);
    atoms_.resize(atoi(line.data()));

    index = 0;
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) { return false; }
        if(!atoms_[i].parse(line)) { return false; }
        
        atomType = searchForIdenticalAtomType(atoms_[i], i-1);
        if(atomType == -1)
        {
            index++;
            atoms_[i].atomType_ = index;
            atoms_[i].firstAtomOfThisType_ = true;
        }
        else
        {
            atoms_[i].atomType_ = atomType;
            atoms_[i].firstAtomOfThisType_ = false;
        }

        atoms_[i].atomTypeGlobal_ = atoms_[i].atomType_ + offsetAtomTypeGlobal;
    }


    // Retrieve comments
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        
    } while((line[0] == '#') && !lastLineInFile);
    
    CASCIIUtility::removeWhiteSpace(line);
    numComments_ = atoi(line.data());
    
    for(int i=0; i<numComments_; i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) { return false; }
    }


    // Retrieve bonds
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        
    } while((line[0] == '#') && !lastLineInFile);
    
    CASCIIUtility::removeWhiteSpace(line);
    bonds_.resize(atoi(line.data()));
    
    for(int i=0; i<(int)bonds_.size(); i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) { return false; }
        if(!bonds_[i].parse(line)) { return false; }
    }
    
    
    // Retrieve angles
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        
    } while((line[0] == '#') && !lastLineInFile);
    
    CASCIIUtility::removeWhiteSpace(line);
    angles_.resize(atoi(line.data()));
    
    for(int i=0; i<(int)angles_.size(); i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) { return false; }
        if(!angles_[i].parse(line)) { return false; }
    }

    
    // Retrieve dihedrals
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        
    } while((line[0] == '#') && !lastLineInFile);
    
    CASCIIUtility::removeWhiteSpace(line);
    dihedrals_.resize(atoi(line.data()));
    
    for(int i=0; i<(int)dihedrals_.size(); i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) { return false; }
        if(!dihedrals_[i].parse(line, -1)) { return false; }
    }
    
    
    // Retrieve optionals
    do
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        CASCIIUtility::removeWhiteSpace(line);
        if(line == "tors1")
        {
            if(!handleTors(fileHandle, 1)) { return false; }
        }
        if(line == "tors5")
        {
            if(!handleTors(fileHandle, 5)) { return false; }
        }
        
    } while(!lastLineInFile);

    // Is this MMOL an fSPC type MMOL? Start from begining of file to find out
    fseek(fileHandle, 0, SEEK_SET);
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        std::string szText = line;
        if(szText.find("fSPC") != std::string::npos)
        {
            flexibleSPC_ = true;
        }
        
    } while(!lastLineInFile);
    
    return true;
}

int CMolDB::CMolecule::searchForIdenticalAtomType(const CAtom &atom, int lastIndex) const
{
    if(lastIndex == -1) lastIndex = (int)atoms_.size() - 1;
    
    for(int i=0; i<(lastIndex+1); i++)
    {
        if(atom == atoms_[i]) return atoms_[i].atomType_;
    }
    
    return -1;
}

int CMolDB::CMolecule::getNumDistinctAtomTypes() const
{
    int count = 0;
    
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        if(atoms_[i].firstAtomOfThisType_) count++;
    }
    
    return count;
}

bool CMolDB::CMolecule::handleTors(FILE* fileHandle, int torsType)
{
    std::string line;
    bool lastLineInFile;
    int lastCount;
     
    do
    {
        line = CFileUtility::readLine(fileHandle, lastLineInFile);
        
    } while((line[0] == '#') && !lastLineInFile);
    
    CASCIIUtility::removeWhiteSpace(line);
    lastCount = (int)dihedrals_.size();
    dihedrals_.resize(lastCount + atoi(line.data()));
    
    for(int i=lastCount; i<(int)dihedrals_.size(); i++)
    {
        do
        {
            line = CFileUtility::readLine(fileHandle, lastLineInFile);
            
        } while((line[0] == '#') && !lastLineInFile);
        
        if(lastLineInFile) return false;
        if(!dihedrals_[i].parse(line, torsType)) return false;
    }
    
    return true;
}


void CMolDB::CMolecule::print() const
{
    printf("Number of atoms in molecule: %i\r\n-----------------------------------------\r\n", (int)atoms_.size());
    for(int i=0; i<(int)atoms_.size(); i++)
    {
        atoms_[i].print();
    }
}


////////////////////////////////////////////////////////////////////////
// CMolDB
////////////////////////////////////////////////////////////////////////

CMolDB::CMolDB()
{
}

bool CMolDB::addMoleculeType(const char* fileName)
{
    FILE* handle = fopen(fileName, "r");
    int offsetAtomTypeGlobal = 0;
    
    if(!handle) return false;
    
    for(int i=0; i<(int)moleculeTypes_.size(); i++)
    {
        offsetAtomTypeGlobal+= moleculeTypes_[i].getNumDistinctAtomTypes();
    }
    
    moleculeTypes_.emplace_back(CMolecule());
    if(!moleculeTypes_[moleculeTypes_.size()-1].load(handle, offsetAtomTypeGlobal))
    {
        fclose(handle);
        return false;
    }
    
    fclose(handle);
    
    return true;
}

void CMolDB::print() const
{
    for(int i=0; i<(int)moleculeTypes_.size(); i++)
    {
        moleculeTypes_[i].print();
        printf("\r\n");
    }
}

