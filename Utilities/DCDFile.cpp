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

#include <string.h>
#include <float.h>
#include "DCDFile.h"

BEGIN_CUDA_COMPATIBLE()

////////////////////////////////////////////////////////////////////////
// CDCDFile::CMainHeader
////////////////////////////////////////////////////////////////////////

CDCDFile::CMainHeader::CMainHeader()
{
    ID_ = "    ";
    descriptA_.resize(80, ' ');
    descriptB_.resize(80, ' ');

    nSets_ = 0;
    initStep_ = 0;
    wrtFreq_ = 0;
    timeStep_ = 0.0;
    nAtoms_ = 0;
}

bool CDCDFile::CMainHeader::read(FILE* handle, int& bytesRead)
{
    int readCnt;
    uint32_t temp;
    
    bytesRead = 0;
    if(!handle) return false;
    
    // Read Fortran record length begin
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read filetype ID
    char ID[5];
    readCnt = (int)fread(ID, 1, 4, handle);
    ID[4] = '\0';
    ID_ = ID;
    bytesRead+= readCnt;
    if(readCnt != 4) { return false; }
    
    // Read number of datasets (i.e. records)
    readCnt = (int)fread(&nSets_, 1, sizeof(nSets_), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(nSets_)) { return false; }
    
    // Read initial step
    readCnt = (int)fread(&initStep_, 1, sizeof(initStep_), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(initStep_)) { return false; }
    
    // Read write frequency of datasets
    readCnt = (int)fread(&wrtFreq_, 1, sizeof(wrtFreq_), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(wrtFreq_)) { return false; }
    
    // Read six '0'
    for(int i=0; i<6; i++)
    {
        readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
        bytesRead+= readCnt;
        if(readCnt != sizeof(temp)) { return false; }
    }
    
    // Read timestep in seconds
    readCnt = (int)fread(&timeStep_, 1, sizeof(timeStep_), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(timeStep_)) { return false; }
    
    // Read one '1', eight '0', one '1', Fortran record length end, Fortran record length begin, read one '2'
    for(int i=0; i<13; i++)
    {
        readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
        bytesRead+= readCnt;
        if(readCnt != sizeof(temp)) { return false; }
    }
    
    // Read description A
    char descript[81];
    readCnt = (int)fread(descript, 1, 80, handle);
    descript[80] = '\0';
    descriptA_ = descript;
    bytesRead+= readCnt;
    if(readCnt != 80) { return false; }
    
    // Read description B
    readCnt = (int)fread(descript, 1, 80, handle);
    descript[80] = '\0';
    descriptB_ = descript;
    bytesRead+= readCnt;
    if(readCnt != 80) { return false; }
    
    // Read Fortran record length end, Fortran record length begin
    for(int i=0; i<2; i++)
    {
        readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
        bytesRead+= readCnt;
        if(readCnt != sizeof(temp)) { return false; }
    }
    
    // Read number of atoms
    readCnt = (int)fread(&nAtoms_, 1, sizeof(nAtoms_), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(nAtoms_)) { return false; }
    
    // Read Fortran record length end
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    bytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    return true;
}

bool CDCDFile::CMainHeader::write(FILE* handle, int& bytesWritten) const
{
    int writeCnt;
    uint32_t temp;
    uint32_t writeSeq1[13] = {1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 84, 164, 2 };
    uint32_t writeSeq2[2] = { 164, 4 };
     
    bytesWritten = 0;
    if(!handle) return false;
    
    // Write Fortran record length begin
    temp = 84;
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write filetype ID
    std::string ID = ID_;
    if(ID.size() != 4) ID.resize(4);
    writeCnt = (int)fwrite(ID.data(), 1, 4, handle);
    bytesWritten+= writeCnt;
    if(writeCnt != 4) { return false; }
    
    // Write number of datasets (i.e. records)
    writeCnt = (int)fwrite(&nSets_, 1, sizeof(nSets_), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(nSets_)) { return false; }
    
    // Write initial step
    writeCnt = (int)fwrite(&initStep_, 1, sizeof(initStep_), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(initStep_)) { return false; }
    
    // Write write frequency of datasets
    writeCnt = (int)fwrite(&wrtFreq_, 1, sizeof(wrtFreq_), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(wrtFreq_)) { return false; }
    
    // Write six '0'
    for(int i=0; i<6; i++)
    {
        temp = 0;
        writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
        bytesWritten+= writeCnt;
        if(writeCnt != sizeof(temp)) { return false; }
    }
    
    // Write timestep in seconds
    writeCnt = (int)fwrite(&timeStep_, 1, sizeof(timeStep_), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(timeStep_)) { return false; }
    
    // Write one '1', eight '0', one '1', Fortran record length end, Fortran record length begin, read one '2'
    for(int i=0; i<13; i++)
    {
        writeCnt = (int)fwrite(&writeSeq1[i], 1, sizeof(writeSeq1[i]), handle);
        bytesWritten+= writeCnt;
        if(writeCnt != sizeof(writeSeq1[i])) { return false; }
    }
    
    // Write description A
    std::string descriptA = descriptA_;
    if(descriptA.size() != 80) descriptA.resize(80);
    writeCnt = (int)fwrite(descriptA.data(), 1, 80, handle);
    bytesWritten+= writeCnt;
    if(writeCnt != 80) { return false; }
    
    // Write description B
    std::string descriptB = descriptB_;
    if(descriptB.size() != 80) descriptB.resize(80);
    writeCnt = (int)fwrite(descriptB.data(), 1, 80, handle);
    bytesWritten+= writeCnt;
    if(writeCnt != 80) { return false; }
    
    // Write Fortran record length end, Fortran record length begin
    for(int i=0; i<2; i++)
    {
        writeCnt = (int)fwrite(&writeSeq2[i], 1, sizeof(writeSeq2[i]), handle);
        bytesWritten+= writeCnt;
        if(writeCnt != sizeof(writeSeq2[i])) { return false; }
    }
    
    // Write number of atoms
    writeCnt = (int)fwrite(&nAtoms_, 1, sizeof(nAtoms_), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(nAtoms_)) { return false; }
    
    // Write Fortran record length end
    temp = 4;
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    bytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    return true;
}


////////////////////////////////////////////////////////////////////////
// CDCDFile::CRecordHeader
////////////////////////////////////////////////////////////////////////

CDCDFile::CRecordHeader::CRecordHeader()
{
    boxX_ = 0.0;
    boxY_ = 0.0;
    boxZ_ = 0.0;
    angleX_ = 0.0;
    angleY_ = 0.0;
    angleZ_ = 0.0;
}

bool CDCDFile::CRecordHeader::read(FILE* handle, int& numBytesRead)
{
    int readCnt;
    uint32_t temp;
    unsigned char boxSize[48];
     
    numBytesRead = 0;
    
    // Read Fortran record length begin
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read 48 bytes containing six eight byte doubles
    readCnt = (int)fread(boxSize, 1, 48, handle);
    numBytesRead+= readCnt;
    if(readCnt != 48) { return false; }
    convertBoxSizeToDoubles(boxSize, boxX_, angleX_, boxY_, angleY_, boxZ_, angleZ_);
    
    // Read Fortran record length end
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    return true;
}

bool CDCDFile::CRecordHeader::write(FILE* handle, int& numBytesWritten) const
{
    int writeCnt;
    uint32_t temp;
    unsigned char boxSize[48];
     
    numBytesWritten = 0;
    
    // Read Fortran record length begin
    temp = 48;
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Read 48 bytes containing six eight byte doubles
    convertDoublesToBoxSize(boxSize, boxX_, angleX_, boxY_, angleY_, boxZ_, angleZ_);
    writeCnt = (int)fwrite(boxSize, 1, 48, handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != 48) { return false; }
    
    // Read Fortran record length end
    temp = 48;
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    return true;
}

void CDCDFile::CRecordHeader::print() const
{
    printf("BoxX: %g\t\tAngleX: %g\r\n", boxX_, angleX_);
    printf("BoxY: %g\t\tAngleY: %g\r\n", boxY_, angleY_);
    printf("BoxZ: %g\t\tAngleZ: %g\r\n", boxZ_, angleZ_);
}

void CDCDFile::CRecordHeader::convertBoxSizeToDoubles(const unsigned char* boxSize, double& boxX, double& angX, double& boxY, double& angY, double& boxZ, double& angZ) const
{
    boxX = ((double*)(boxSize))[0];
    angX = ((double*)(boxSize))[1];
    boxY = ((double*)(boxSize))[2];
    angY = ((double*)(boxSize))[3];
    angZ = ((double*)(boxSize))[4];
    boxZ = ((double*)(boxSize))[5];
}

void CDCDFile::CRecordHeader::convertDoublesToBoxSize(unsigned char* boxSize, const double &boxX, const double &angX, const double &boxY, const double &angY, const double &boxZ, const double &angZ) const
{
    ((double*)(boxSize))[0] = boxX;
    ((double*)(boxSize))[1] = angX;
    ((double*)(boxSize))[2] = boxY;
    ((double*)(boxSize))[3] = angY;
    ((double*)(boxSize))[4] = angZ;
    ((double*)(boxSize))[5] = boxZ;
}


////////////////////////////////////////////////////////////////////////
// CDCDFile::CRecord
////////////////////////////////////////////////////////////////////////

CDCDFile::CRecord::CRecord()
{
    numCoordinates_ = 0;
}

CDCDFile::CRecord::~CRecord()
{
}

bool CDCDFile::CRecord::read(FILE* handle, int& numBytesRead)
{
    int readCnt;
    uint32_t temp;
    
    if(!recordHeader_.read(handle, numBytesRead)) return false;
    
    // Read Fortran record length begin
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read X-positions
    xPositions_.resize(numCoordinates_);
    readCnt = (int)fread(xPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesRead+= readCnt * sizeof(float);
    if(readCnt != numCoordinates_) { return false; }
    
    // Read Fortran record length end
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read Fortran record length begin
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read Y-positions
    yPositions_.resize(numCoordinates_);
    readCnt = (int)fread(yPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesRead+= readCnt * sizeof(float);
    if(readCnt != numCoordinates_) { return false; }
    
    // Read Fortran record length end
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read Fortran record length begin
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }
    
    // Read Z-positions
    zPositions_.resize(numCoordinates_);
    readCnt = (int)fread(zPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesRead+= readCnt * sizeof(float);
    if(readCnt != numCoordinates_) { return false; }
    
    // Read Fortran record length end
    readCnt = (int)fread(&temp, 1, sizeof(temp), handle);
    numBytesRead+= readCnt;
    if(readCnt != sizeof(temp)) { return false; }

    return true;
}

bool CDCDFile::CRecord::write(FILE* handle, int& numBytesWritten) const
{
    int writeCnt;
    uint32_t temp;
    
    if(!recordHeader_.write(handle, numBytesWritten)) return false;

    // Write Fortran record length begin
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write X-positions
    if(xPositions_.size() != numCoordinates_) return false;
    writeCnt = (int)fwrite(xPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesWritten+= writeCnt * sizeof(float);
    if(writeCnt != numCoordinates_) { return false; }
    
    // Write Fortran record length end
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write Fortran record length begin
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write Y-positions
    if(yPositions_.size() != numCoordinates_) return false;
    writeCnt = (int)fwrite(yPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesWritten+= writeCnt * sizeof(float);
    if(writeCnt != numCoordinates_) { return false; }
    
    // Write Fortran record length end
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write Fortran record length begin
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }
    
    // Write Z-positions
    if(zPositions_.size() != numCoordinates_) return false;
    writeCnt = (int)fwrite(zPositions_.data(), sizeof(float), numCoordinates_, handle);
    numBytesWritten+= writeCnt * sizeof(float);
    if(writeCnt != numCoordinates_) { return false; }
    
    // Write Fortran record length end
    temp = numCoordinates_ * sizeof(float);
    writeCnt = (int)fwrite(&temp, 1, sizeof(temp), handle);
    numBytesWritten+= writeCnt;
    if(writeCnt != sizeof(temp)) { return false; }

    return true;
}

void CDCDFile::CRecord::init(uint32_t numCoordinates, bool resizeArray)
{
    numCoordinates_ = numCoordinates;
    if(resizeArray)
    {
        xPositions_.resize(numCoordinates);
        yPositions_.resize(numCoordinates);
        zPositions_.resize(numCoordinates);
    }
}

void CDCDFile::CRecord::setPos(uint32_t coordIndex, double x, double y, double z)
{        
    if(coordIndex < xPositions_.size()) xPositions_[coordIndex] = (float)x;
    if(coordIndex < yPositions_.size()) yPositions_[coordIndex] = (float)y;
    if(coordIndex < zPositions_.size()) zPositions_[coordIndex] = (float)z;
}

C3DVector CDCDFile::CRecord::getPos(uint32_t coordIndex) const
{
    C3DVector v;

    size_t sizeX = xPositions_.size();
    size_t sizeY = yPositions_.size();
    size_t sizeZ = zPositions_.size();
    if((coordIndex >= sizeX) || (coordIndex >= sizeY) || (coordIndex >= sizeZ))
    {
        printf("Warning: DCD coordinate index %i is outside one of the number of indices (%i, %i, %i)!\r\n", coordIndex, (int)sizeX, (int)sizeY, (int)sizeZ);
        return v;
    }
    
    v.x_ = (double)xPositions_[coordIndex];
    v.y_ = (double)yPositions_[coordIndex];
    v.z_ = (double)zPositions_[coordIndex];
    
    return v;
}

void CDCDFile::CRecord::print() const
{
    recordHeader_.print();

    size_t sizeX = xPositions_.size();
    size_t sizeY = yPositions_.size();
    size_t sizeZ = zPositions_.size();
    if((numCoordinates_ != sizeX) || (numCoordinates_ != sizeY) || (numCoordinates_ != sizeZ))
    {
        printf("Warning: One of DCD coordinate indices (%i, %i, %i) is unequal to the number of coordinates %i found in file header!\r\n", (int)sizeX, (int)sizeY, (int)sizeZ, numCoordinates_);
        return;
    }

    printf("\r\nNum. Coordinates: %i\r\n--------------------------------\r\n", numCoordinates_);
    for(uint32_t i=0; i<numCoordinates_; i++)
    {
        printf("(%g, %g, %g)\r\n", (double)xPositions_[i], (double)yPositions_[i], (double)zPositions_[i]);
    }
    printf("\r\n");
}


////////////////////////////////////////////////////////////////////////
// CDCDFile
////////////////////////////////////////////////////////////////////////

CDCDFile::CDCDFile()
{
    file_ = nullptr;
    stride_ = 1;
    ptToFirstRecord_ = 0;
}

CDCDFile::~CDCDFile()
{
    close();
}

bool CDCDFile::open(std::string fileName, int stride)
{
    file_ = fopen(fileName.data(), "r");
    if(!file_)
    {
        printf("Error opening file %s!\r\n", fileName.data());
        return false;
    }
    
    stride_ = stride;
    
    int bytesReadHeader;
    if(!mainHeader_.read(file_, bytesReadHeader))
    {
        printf("Error reading DCD header!\r\n");
        fclose(file_);
        return false;
    }
    ptToFirstRecord_ = ftell(file_);
    
    currentRecord_.init(mainHeader_.nAtoms_);
    
    return true;
}

void CDCDFile::saveCurrentRecord(FILE* destFile) const
{
    int numBytesWritten;
    currentRecord_.write(destFile, numBytesWritten);
}

void CDCDFile::close()
{
    if(file_) fclose(file_);

    file_ = nullptr;
    stride_ = 1;
}

void CDCDFile::gotoRecord(int recordIndex)
{
    int numBytesRead;
    
    long pos = ptToFirstRecord_ + (80 + 3*sizeof(float)*mainHeader_.nAtoms_)*recordIndex;
    
    if(!file_)
    {
        printf("Error jumping to DCD record number %i!\r\n", recordIndex);
        return;
    }
    
    fseek(file_, pos, SEEK_SET);
    currentRecord_.read(file_, numBytesRead);
}

void CDCDFile::gotoNextRecord(bool&)
{
    int numBytesRead;
    
    long pos = ptToFirstRecord_ + (80 + 3*sizeof(float)*mainHeader_.nAtoms_);

    if(!file_)
    {
        printf("Error jumping to next DCD record number!\r\n");
        return;
    }

    fseek(file_, pos, SEEK_SET);
    currentRecord_.read(file_, numBytesRead);
}

C3DRect CDCDFile::getCurrentBoundingBox() const
{
    C3DRect boundingBox;

    boundingBox.rLow_.x_ = DBL_MAX;
    boundingBox.rLow_.y_ = DBL_MAX;
    boundingBox.rLow_.z_ = DBL_MAX;
    
    boundingBox.rHigh_.x_ = -DBL_MAX;
    boundingBox.rHigh_.y_ = -DBL_MAX;
    boundingBox.rHigh_.z_ = -DBL_MAX;

    unsigned long numCoordinates = currentRecord_.getNumCoordinates();
    std::pair<const float*, size_t> xPositions = currentRecord_.getCoordinateDataX();
    std::pair<const float*, size_t> yPositions = currentRecord_.getCoordinateDataY();
    std::pair<const float*, size_t> zPositions = currentRecord_.getCoordinateDataZ();
    for(unsigned long i=0; i<numCoordinates; i++)
    {
        if(i >= xPositions.second) continue;
        if(i >= yPositions.second) continue;
        if(i >= zPositions.second) continue;

        if((double)xPositions.first[i] < boundingBox.rLow_.x_) boundingBox.rLow_.x_ = (double)xPositions.first[i];
        if((double)yPositions.first[i] < boundingBox.rLow_.y_) boundingBox.rLow_.y_ = (double)yPositions.first[i];
        if((double)zPositions.first[i] < boundingBox.rLow_.z_) boundingBox.rLow_.z_ = (double)zPositions.first[i];

        if((double)xPositions.first[i] > boundingBox.rHigh_.x_) boundingBox.rHigh_.x_ = (double)xPositions.first[i];
        if((double)yPositions.first[i] > boundingBox.rHigh_.y_) boundingBox.rHigh_.y_ = (double)yPositions.first[i];
        if((double)zPositions.first[i] > boundingBox.rHigh_.z_) boundingBox.rHigh_.z_ = (double)zPositions.first[i];
    }
    
    return boundingBox;
}

C3DRect CDCDFile::getCurrentPBC() const
{
    C3DRect boundingBox = getCurrentBoundingBox();
    C3DVector center = boundingBox.getCenter();
    C3DRect pbc;
    CRecordHeader recordHeader = currentRecord_.getRecordHeader();
    double widthX_2 = recordHeader.boxX_ / 2.0;
    double widthY_2 = recordHeader.boxY_ / 2.0;
    double widthZ_2 = recordHeader.boxZ_ / 2.0;

    pbc.rLow_.x_ = center.x_ - widthX_2;
    pbc.rHigh_.x_ = center.x_ + widthX_2;

    pbc.rLow_.y_ = center.y_ - widthY_2;
    pbc.rHigh_.y_ = center.y_ + widthY_2;

    pbc.rLow_.z_ = center.z_ - widthZ_2;
    pbc.rHigh_.z_ = center.z_ + widthZ_2;
    
    return pbc;
}

void CDCDFile::setCoordinate(int coordinateIndex, int coordinate, double val)
{
    std::pair<const float*, size_t> xPositions = currentRecord_.getCoordinateDataX();
    std::pair<const float*, size_t> yPositions = currentRecord_.getCoordinateDataY();
    std::pair<const float*, size_t> zPositions = currentRecord_.getCoordinateDataZ();

    if((coordinate == 0) && (coordinate < (int)xPositions.second))  currentRecord_.setXPos(coordinateIndex, val);
    if((coordinate == 1) && (coordinate < (int)yPositions.second))  currentRecord_.setYPos(coordinateIndex, val);
    if((coordinate == 2) && (coordinate < (int)zPositions.second))  currentRecord_.setZPos(coordinateIndex, val);
}

END_CUDA_COMPATIBLE()
