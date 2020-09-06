#pragma once
#include <stdio.h>
#include <string>
#include <vector>
#include "3DVector.h"
#include "3DRect.h"

BEGIN_CUDA_COMPATIBLE()

#ifndef _UINT32_T
#define _UINT32_T
typedef unsigned int uint32_t;
#endif

class CDCDFile
{
public:
    struct CMainHeader
    {
    public:
        CMainHeader();
        
    public:
        bool read(FILE* handle, int& bytesRead);
        bool write(FILE* handle, int& bytesWritten) const;
        
    public:
        std::string ID_;
        uint32_t nSets_;
        uint32_t initStep_;
        uint32_t wrtFreq_;
        float timeStep_;
        std::string descriptA_;
        std::string descriptB_;
        uint32_t nAtoms_;
    };
    
    struct CRecordHeader
    {
    public:
        CRecordHeader();
        
    public:
        bool read(FILE* handle, int& numBytesRead);
        bool write(FILE* handle, int& numBytesWritten) const;
        void print() const;
        
    private:
        void convertBoxSizeToDoubles(const unsigned char* boxSize, double& boxX, double& angX, double& boxY, double& angY, double& boxZ, double& angZ) const;
        void convertDoublesToBoxSize(unsigned char* boxSize, const double& boxX, const double& angX, const double& boxY, const double& angY, const double& boxZ, const double& angZ) const;
        
    public:
        double boxX_;
        double boxY_;
        double boxZ_;
        double angleX_;
        double angleY_;
        double angleZ_;
    };
    
    class CRecord
    {
    public:
        CRecord();
        ~CRecord();
        
    public:
        bool read(FILE* handle, int& numBytesRead);
        bool write(FILE* handle, int& numBytesWritten) const;
        void init(uint32_t numCoordinates, bool resizeArray=false);
        void setPos(uint32_t coordIndex, double x, double y, double z);
        void setXPos(uint32_t coordIndex, double x) { xPositions_[coordIndex] = (float)x; }
        void setYPos(uint32_t coordIndex, double y) { yPositions_[coordIndex] = (float)y; }
        void setZPos(uint32_t coordIndex, double z) { zPositions_[coordIndex] = (float)z; }
        std::pair<const float*, size_t> getCoordinateDataX() const { return std::pair<const float*, size_t>(xPositions_.data(), xPositions_.size()); }
        std::pair<const float*, size_t> getCoordinateDataY() const { return std::pair<const float*, size_t>(yPositions_.data(), yPositions_.size()); }
        std::pair<const float*, size_t> getCoordinateDataZ() const { return std::pair<const float*, size_t>(zPositions_.data(), zPositions_.size()); }
        C3DVector getPos(uint32_t coordIndex) const;
        void setRecordHeader(const CRecordHeader& recordHeader) { recordHeader_ = recordHeader; }
        CRecordHeader getRecordHeader() const { return recordHeader_; }
        void print() const;
        uint32_t getNumCoordinates() const { return numCoordinates_; }
        
    private:
        CRecordHeader recordHeader_;
        std::vector<float> xPositions_;
        std::vector<float> yPositions_;
        std::vector<float> zPositions_;

    private:
        uint32_t numCoordinates_;
    };
    
public:
    CDCDFile();
    ~CDCDFile();
    
public:
    bool open(std::string fileName, int stride=1);
    void saveCurrentRecord(FILE* destFile) const;
    void close();
    void gotoRecord(int recordIndex);
    void gotoNextRecord(bool& lastRecordInFile);
    int getNumRecords() const { return mainHeader_.nSets_; }
    int getNumCoordinatesInRecord() const { return (int)currentRecord_.getNumCoordinates(); }
    C3DVector getCoordinate(int coordinateIndex) const { return currentRecord_.getPos(coordinateIndex); }
    void setCoordinate(int coordinateIndex, C3DVector r) { currentRecord_.setPos(coordinateIndex, r.x_, r.y_, r.z_); }
    void setCoordinate(int coordinateIndex, int coordinate, double val);
    void setXCoordinate(int coordinateIndex, double x) { currentRecord_.setXPos((uint32_t)coordinateIndex, x); }
    void setYCoordinate(int coordinateIndex, double y) { currentRecord_.setYPos((uint32_t)coordinateIndex, y); }
    void setZCoordinate(int coordinateIndex, double z) { currentRecord_.setZPos((uint32_t)coordinateIndex, z); }
    const float* getCoordinateDataX() { return currentRecord_.getCoordinateDataX().first; }
    const float* getCoordinateDataY() { return currentRecord_.getCoordinateDataY().first; }
    const float* getCoordinateDataZ() { return currentRecord_.getCoordinateDataZ().first; }
    CMainHeader getDCDHeader() const { return mainHeader_; }
    CRecordHeader getCurrentRecordHeader() const { return currentRecord_.getRecordHeader(); }
    C3DRect getCurrentBoundingBox() const;
    C3DRect getCurrentPBC() const;
    
private:
    FILE* file_;
    int stride_;
    CMainHeader mainHeader_;
    CRecord currentRecord_;
    long ptToFirstRecord_;
};

END_CUDA_COMPATIBLE()
