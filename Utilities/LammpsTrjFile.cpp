//
// Copyright (C) 2025 Richard Olsen.
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
#include <fstream>
#include "ASCIIUtility.h"
#include "LammpsTrjFile.h"

BEGIN_CUDA_COMPATIBLE()

CLammpsTrjFile::~CLammpsTrjFile()
{
    close();
}

bool CLammpsTrjFile::open(std::string fileName)
{
    // Open File
    file_ = std::make_shared<std::ifstream>(fileName);
    if(!file_ || !file_->is_open())
    {
        printf("Error opening file %s!\r\n", fileName.data());
        return false;
    }

    // Map all records to file positions
    const std::streampos startPos = file_->tellg();

    bool successfulLine = false;
    std::string line;
    do
    {
        const std::streampos filePos = file_->tellg();
        successfulLine = std::getline(*file_, line) ? true : false;
        if(successfulLine && (CASCIIUtility::findString("TIMESTEP", line) > -1))
        {
            recordFilePositions_.emplace_back(filePos);
        }

    } while(successfulLine);

    file_->clear();
    file_->seekg(startPos);

    // Load the first record
    gotoRecord(0);
    return true;
}

void CLammpsTrjFile::close()
{
    if(file_) file_->close();

    file_ = nullptr;
}

void CLammpsTrjFile::gotoRecord(int recordIndex)
{
    if(recordIndex < (int)recordFilePositions_.size())
    {
        file_->clear();
        file_->seekg(recordFilePositions_[recordIndex]);
        currentRecord_ = std::make_shared<CRecord>(file_);
    }
}

int CLammpsTrjFile::getNumRecords() const
{
    return (int)recordFilePositions_.size();
}

int CLammpsTrjFile::getNumCoordinatesInRecord() const
{
    if(!currentRecord_) return -1;
    return currentRecord_->getNumAtoms();
}

int CLammpsTrjFile::getTimeStepIndex() const
{
    return currentRecord_->getTimeStepIndex();
}

C3DVector CLammpsTrjFile::getCoordinate(int atomIndex) const
{
    if(!currentRecord_) return C3DVector();
    return currentRecord_->getCoordinate(atomIndex);
}

C3DRect CLammpsTrjFile::getCurrentPBC() const
{
    if(!currentRecord_) return C3DRect();
    return currentRecord_->getPBC();
}



CLammpsTrjFile::CAtom::CAtom(const std::string& line, const int& indexAtom, const int& indexType, const int& indexX, const int& indexY, const int& indexZ)
{
    if(!parse(line, indexAtom, indexType, indexX, indexY, indexZ))
    {
        printf("Error: Could not parse atom!\r\n");
    }
}

void CLammpsTrjFile::CAtom::unscaleX(const std::pair<double, double>& bbXSpan)
{
    const double xSpan = bbXSpan.second - bbXSpan.first;
    x_ = bbXSpan.first + xSpan*x_;
}

void CLammpsTrjFile::CAtom::unscaleY(const std::pair<double, double>& bbYSpan)
{
    const double ySpan = bbYSpan.second - bbYSpan.first;
    y_ = bbYSpan.first + ySpan*y_;
}

void CLammpsTrjFile::CAtom::unscaleZ(const std::pair<double, double>& bbZSpan)
{
    const double zSpan = bbZSpan.second - bbZSpan.first;
    z_ = bbZSpan.first + zSpan*z_;
}

int CLammpsTrjFile::CAtom::getAtomIndex() const
{
    return atomNumber_ - 1;
}

C3DVector CLammpsTrjFile::CAtom::getPosition() const
{
    return C3DVector(x_, y_, z_);
}

bool CLammpsTrjFile::CAtom::parse(const std::string& line, const int& indexAtom, const int& indexType, const int& indexX, const int& indexY, const int& indexZ)
{
    std::vector<std::string> words = CASCIIUtility::getWords(line);
    const int numWords = (int)words.size();
    if((indexAtom >= 0) && (indexAtom >= numWords)) return false;
    if((indexType >= 0) && (indexType >= numWords)) return false;
    if((indexX >= 0) && (indexX >= numWords)) return false;
    if((indexY >= 0) && (indexY >= numWords)) return false;
    if((indexZ >= 0) && (indexZ >= numWords)) return false;

    atomNumber_ = atoi(words[indexAtom].data());
    typeIndex_ = atoi(words[indexType].data());
    x_ = atof(words[indexX].data());
    y_ = atof(words[indexY].data());
    z_ = atof(words[indexZ].data());

    return true;
}



CLammpsTrjFile::CRecord::CRecord(std::shared_ptr<std::ifstream> file)
{
    if(!read(file))
    {
        printf("Error: could not read file!\r\n");
    }
}

bool CLammpsTrjFile::CRecord::read(std::shared_ptr<std::ifstream> file)
{
    if(!file) return false;

    std::string line;
    while(std::getline(*file, line))
    {
        if((CASCIIUtility::findString("TIMESTEP", line) > -1) && std::getline(*file, line))
        {
            timeStep_ = atoi(line.data());
        }
        else if((CASCIIUtility::findString("NUMBER OF ATOMS", line) > -1) && std::getline(*file, line))
        {
            numAtoms_ = atoi(line.data());
            atoms_.resize(numAtoms_);
        }
        else if(CASCIIUtility::findString("BOX BOUNDS", line) > -1)
        {
            if(CASCIIUtility::findString("xy xz yz", line) > -1)
            {
                printf("Error: MolTwister does not yet support triclinic tilted simulation boxes. Found triclinic box xy xz yz flags in LAMMPS trajectory file. If these are set and system is not triclinic, you may attempt to remove these flags and the corresponding values from the *.lammpstrj file manually!\r\n");
                return false;
            }

            if(std::getline(*file, line))
            {
                double low = atof(CASCIIUtility::getWord(line, 0).data());
                double high = atof(CASCIIUtility::getWord(line, 1).data());
                pbcXSpan_ = std::pair<double, double>(low, high);
            }
            if(std::getline(*file, line))
            {
                double low = atof(CASCIIUtility::getWord(line, 0).data());
                double high = atof(CASCIIUtility::getWord(line, 1).data());
                pbcYSpan_ = std::pair<double, double>(low, high);
            }
            if(std::getline(*file, line))
            {
                double low = atof(CASCIIUtility::getWord(line, 0).data());
                double high = atof(CASCIIUtility::getWord(line, 1).data());
                pbcZSpan_ = std::pair<double, double>(low, high);
            }
        }
        else if(CASCIIUtility::findString("ATOMS", line) > -1)
        {
            int indexAtom = -1;
            int indexType = -1;
            int indexX = -1;
            int indexY = -1;
            int indexZ = -1;
            int indexXS = -1;
            int indexYS = -1;
            int indexZS = -1;
            int indexXU = -1;
            int indexYU = -1;
            int indexZU = -1;

            std::vector<std::string> words = CASCIIUtility::getWords(line);
            for(int i=0; i<(int)words.size(); i++)
            {
                words[i] = CASCIIUtility::trimFromLeft(words[i]);
                words[i] = CASCIIUtility::trimFromRight(words[i]);

                if(words[i] == "id") indexAtom = i - 2;
                if(words[i] == "type") indexType = i - 2;
                if(words[i] == "x") indexX = i - 2;
                if(words[i] == "y") indexY = i - 2;
                if(words[i] == "z") indexZ = i - 2;
                if(words[i] == "xs") indexXS = i - 2;
                if(words[i] == "ys") indexYS = i - 2;
                if(words[i] == "zs") indexZS = i - 2;
                if(words[i] == "xu") indexXU = i - 2;
                if(words[i] == "yu") indexYU = i - 2;
                if(words[i] == "zu") indexZU = i - 2;
            }

            bool unscaleX = false;
            if(indexX > -1) { } // keep indexX
            else if (indexXS > -1) { indexX = indexXS; unscaleX = true; }
            else if (indexXU > -1) { indexX = indexXU; }
            else return false;

            bool unscaleY = false;
            if(indexY > -1) { } // keep indexY
            else if (indexYS > -1) { indexY = indexYS; unscaleY = true; }
            else if (indexYU > -1) { indexY = indexYU; }
            else return false;

            bool unscaleZ = false;
            if(indexZ > -1) { } // keep indexZ
            else if (indexZS > -1) { indexZ = indexZS; unscaleZ = true; }
            else if (indexZU > -1) { indexZ = indexZU; }
            else return false;

            for(int i=0; i<numAtoms_; i++)
            {
                if(std::getline(*file, line))
                {
                    CAtom atom(line, indexAtom, indexType, indexX, indexY, indexZ);
                    if(unscaleX) atom.unscaleX(pbcXSpan_);
                    if(unscaleY) atom.unscaleY(pbcYSpan_);
                    if(unscaleZ) atom.unscaleZ(pbcZSpan_);

                    const int indexOfAtom = atom.getAtomIndex();
                    if(indexOfAtom < numAtoms_)
                    {
                        atoms_[indexOfAtom] = atom;
                    }
                }
            }

            break;
        }
    }

    return true;
}

int CLammpsTrjFile::CRecord::getNumAtoms() const
{
    return numAtoms_;
}

int CLammpsTrjFile::CRecord::getTimeStepIndex() const
{
    return timeStep_;
}

C3DVector CLammpsTrjFile::CRecord::getCoordinate(const int& index) const
{
    return atoms_[index].getPosition();
}

C3DRect CLammpsTrjFile::CRecord::getPBC() const
{
    return C3DRect(C3DVector(pbcXSpan_.first, pbcYSpan_.first, pbcZSpan_.first),
                   C3DVector(pbcXSpan_.second, pbcYSpan_.second, pbcZSpan_.second));
}

END_CUDA_COMPATIBLE()
