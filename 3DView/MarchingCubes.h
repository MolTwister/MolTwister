#pragma once
////////////////////////////////////////////////////////////////////////////////////////////////
// Code based on an implementation by Paul Bourke ( http://paulbourke.net/geometry/polygonise/ )
////////////////////////////////////////////////////////////////////////////////////////////////

#include "Utilities/3DVector.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMarchingCubesVoxel
{
public:
    CMarchingCubesVoxel();
    CMarchingCubesVoxel(C3DVector* R, double* val);
    
public:
    C3DVector R_[8];
    C3DVector N_[8];
    double val_[8];
};

class CMarchingCubesTriangle
{
public:
    CMarchingCubesTriangle();
    CMarchingCubesTriangle(C3DVector* points);
    
public:
    C3DVector R_[3];
    C3DVector N_[3];
};

class CMarchingCubes
{
public:
    CMarchingCubes();
    
public:
    static int calcTriangles(const CMarchingCubesVoxel& voxel, double isoLevel, CMarchingCubesTriangle* triangles);
    static int calcTrianglesSmooth(const CMarchingCubesVoxel& voxel, double isoLevel, CMarchingCubesTriangle* triangles);

private:
    static C3DVector vertexInterp(double isoLevel, const C3DVector& p1, const C3DVector& p2, const double& val1, const double& val2);
    
private:
    static int edgeTable_[256];
    static int triTable_[256][16];
};

END_CUDA_COMPATIBLE()
