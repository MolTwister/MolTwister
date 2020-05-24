#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "../../../CudaDefinitions.h"
#include "../../../Cmd/Tools/CudaDeviceList.h"
#include "../../../Cmd/MolTwisterCmdMolDyn/Integrators/Particle3D.h"
#include "CudaTests_cuda.h"

#define ARRAY_SIZE 10
#define SIZE_M 5
#define SIZE_N 6

class CInfo_cuda
{
public:
    HOSTDEV_CALLABLE CInfo_cuda() { info1_ = 0; info2_ = 0.0f; }
    HOSTDEV_CALLABLE CInfo_cuda(int info1, float info2) { info1_ = info1; info2_ = info2; }

public:
    HOSTDEV_CALLABLE int getInfo1() { return info1_; }
    HOSTDEV_CALLABLE float getInfo2() { return info2_; }

private:
    int info1_;
    float info2_;
};

class CParticle3D_cuda
{
public:
    HOSTDEV_CALLABLE CParticle3D_cuda() { i_ = -1; m_ = x_[0] = x_[1] = x_[2] = p_[0] = p_[1] = p_[2] = 0.0f; }
    HOSTDEV_CALLABLE CParticle3D_cuda(int i, float m, float x, float y, float z, float px, float py, float pz)
    {
        i_ = i;

        m_ = m;

        x_[0] = x;
        x_[1] = y;
        x_[2] = z;

        p_[0] = px;
        p_[1] = py;
        p_[2] = pz;
    }

public:
    int i_;
    float m_;
    float x_[3];
    float p_[3];
    CInfo_cuda info_;
};

CUDA_GLOBAL void kernelAddBIntoA(int* A, int* B)
{
    long lIdx = blockDim.x*blockIdx.x + threadIdx.x;
    if(lIdx < ARRAY_SIZE)
    {
        A[lIdx] = A[lIdx] + B[lIdx];
    }
}

HOSTDEV_CALLABLE int toIndex(int row, int column, int columnCount)
{
    return columnCount*row + column;
}

class TransfStep
{
public:
    HOSTDEV_CALLABLE TransfStep(CInfo_cuda* info2DMatrix, int rows, int columns)
    {
        rowsinfo2DMatrix_ = rows;
        columnsinfo2DMatrix_ = columns;
        info2DMatrix_ = info2DMatrix;
    }

public:
    HOSTDEV_CALLABLE CParticle3D_cuda operator()(const CParticle3D_cuda& in)
    {
        CParticle3D_cuda out;

        out.m_ = in.m_ / 2.0f;
        for(size_t i=0; i<3; i++)
        {
            out.x_[i] = in.x_[i] / 4.0f;
            out.p_[i] = in.p_[i] / 8.0f;
        }

        if(in.i_ < (rowsinfo2DMatrix_*columnsinfo2DMatrix_)) out.info_ = info2DMatrix_[in.i_];

        return out;
    }

private:
    CInfo_cuda* info2DMatrix_;
    int rowsinfo2DMatrix_;
    int columnsinfo2DMatrix_;
};

void CCudaTest_cuda::addBIntoA(int* A, int* B)
{
    std::pair<long, long> NBlocks = CCudaDeviceList().CalcNumGPUBlocks1D(ARRAY_SIZE);

    // Upload values from A onto the device (GPU)
    int* devA;
    cudaMalloc(&devA, ARRAY_SIZE*sizeof(int));
    cudaMemcpy(devA, A, ARRAY_SIZE*sizeof(int), cudaMemcpyHostToDevice);

    // Upload values from B onto the device (GPU)
    int* devB;
    cudaMalloc(&devB, ARRAY_SIZE*sizeof(int));
    cudaMemcpy(devB, B, ARRAY_SIZE*sizeof(int), cudaMemcpyHostToDevice);

    // Perform parallel manipulations of device versions of A and B vectors
    kernelAddBIntoA<<<NBlocks.first,NBlocks.second>>>(devA, devB);

    // Copy results, stored in device A, from the device to the host
    cudaMemcpy(A, devA, ARRAY_SIZE*sizeof(int), cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    // Free deveice allocated memory
    cudaFree(devA);
    cudaFree(devB);
}

void CCudaTest_cuda::testModifyAtomList(FILE* stdOut)
{
    // Create some common information, organized into a MxN array, to be shared inside the thrust::transform() call
    thrust::host_vector<CInfo_cuda> hostInfo2DMatrix(SIZE_M * SIZE_N);
    int info1 = 2;
    float info2 = 3.5;
    for(size_t r=0; r<SIZE_M; r++)
    {
        for(size_t c=0; c<SIZE_N; c++)
        {
            hostInfo2DMatrix[toIndex((int)r, (int)c, SIZE_N)] = CInfo_cuda(info1, info2);
            info1+= 2;
            info2+= 1.3;
        }
    }

    // Upload common information to device and create a pointer that can be used on the host
    thrust::device_vector<CInfo_cuda> kernelInfo2DMatrix = hostInfo2DMatrix;
    CInfo_cuda* kernelInfo2DMatrixPtr = thrust::raw_pointer_cast(&kernelInfo2DMatrix[0]);

    // Create a vector of particles on the host and set the particle values
    thrust::host_vector<CParticle3D_cuda> hostAtomList(ARRAY_SIZE);

    for(size_t i=0; i<hostAtomList.size(); i++)
    {
        hostAtomList[i] = CParticle3D_cuda(i,
                                           float(i),
                                           float(i)*1.0f, float(i)*2.0f, float(i)*3.0,
                                           float(i)*4.0f, float(i)*5.0f, float(i)*6.0);
    }

    // Create a corresponding vector on the device and copy the host content to the device
    thrust::device_vector<CParticle3D_cuda> kernelAtomList;
    kernelAtomList = hostAtomList;

    // Run a transformation on each particle, on the device (in parallel), defined by the functor, TransfStep
    TransfStep func(kernelInfo2DMatrixPtr, SIZE_M, SIZE_N);
    thrust::device_vector<CParticle3D_cuda> kernelOutput(ARRAY_SIZE);
    thrust::transform(kernelAtomList.begin(), kernelAtomList.end(), kernelOutput.begin(), func);

    // Create a host vector for the result and copy the device content result to it
    thrust::host_vector<CParticle3D_cuda> hostOutput;
    hostOutput = kernelOutput;

    // Print the results
    fprintf(stdOut, "\tPerform a modification (divide by 2, 4 and 8) of mass, m, position, x, and momentum p, of a list of atoms, using CUDA Thrust (info should be (2+i*2, 3.5+i*1.3))\r\n");
    for(size_t i=0; i<ARRAY_SIZE; i++)
    {
        fprintf(stdOut, "\t%i: m_in = %.2f; x_in = (%.2f, %.2f, %.2f); p_in = (%.2f, %.2f, %.2f); (m_in/2, x_in/4, p_in/8) = (%.2f, (%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f)); info=(%i, %.4f)\r\n",
                int(i),
                hostAtomList[i].m_,
                hostAtomList[i].x_[0], hostAtomList[i].x_[1], hostAtomList[i].x_[2],
                hostAtomList[i].p_[0], hostAtomList[i].p_[1], hostAtomList[i].p_[2],
                hostOutput[i].m_,
                hostOutput[i].x_[0], hostOutput[i].x_[1], hostOutput[i].x_[2],
                hostOutput[i].p_[0], hostOutput[i].p_[1], hostOutput[i].p_[2],
                hostOutput[i].info_.getInfo1(), hostOutput[i].info_.getInfo2());
    }
}

int CCudaTest_cuda::getArraySize()
{
    return ARRAY_SIZE;
}
