#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "../../../CudaDefinitions.h"
#include "../../../Cmd/Tools/CudaDeviceList.h"
#include "../../../Cmd/MolTwisterCmdMolDyn/Integrators/Particle3D.h"
#include "CudaTests_cuda.h"

#define ARRAY_SIZE 10

class CParticle3D_cuda
{
public:
    HOSTDEV_CALLABLE CParticle3D_cuda() { m_ = x_[0] = x_[1] = x_[2] = p_[0] = p_[1] = p_[2] = 0.0f; }
    HOSTDEV_CALLABLE CParticle3D_cuda(float m, float x, float y, float z, float px, float py, float pz)
    {
        m_ = m;

        x_[0] = x;
        x_[1] = y;
        x_[2] = z;

        p_[0] = px;
        p_[1] = py;
        p_[2] = pz;
    }

public:
    float m_;
    float x_[3];
    float p_[3];
};

CUDA_GLOBAL void kernelAddBIntoA(int* A, int* B)
{
    long lIdx = blockDim.x*blockIdx.x + threadIdx.x;
    if(lIdx < ARRAY_SIZE)
    {
        A[lIdx] = A[lIdx] + B[lIdx];
    }
}

struct transfStep
{
    HOSTDEV_CALLABLE
    CParticle3D_cuda operator()(const CParticle3D_cuda& in)
    {
        CParticle3D_cuda out;

        out.m_ = in.m_ / 2.0f;
        for(size_t i=0; i<3; i++)
        {
            out.x_[i] = in.x_[i] / 4.0f;
            out.p_[i] = in.p_[i] / 8.0f;
        }

        return out;
    }
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
    // Create a vector of particles on the host and set the particle values
    thrust::host_vector<CParticle3D_cuda> hostAtomList(ARRAY_SIZE);

    for(size_t i=0; i<hostAtomList.size(); i++)
    {
        hostAtomList[i] = CParticle3D_cuda(float(i),
                                           float(i)*1.0f, float(i)*2.0f, float(i)*3.0,
                                           float(i)*4.0f, float(i)*5.0f, float(i)*6.0);
    }

    // Create a corresponding vector on the device and copy the host content to the device
    thrust::device_vector<CParticle3D_cuda> kernelAtomList;
    kernelAtomList = hostAtomList;

    // Run a transformation on each particle, on the device (in parallel), defined by the functor, TransfStep
    transfStep func;
    thrust::device_vector<CParticle3D_cuda> kernelOutput(ARRAY_SIZE);
    thrust::transform(kernelAtomList.begin(), kernelAtomList.end(), kernelOutput.begin(), func);

    // Create a host vector for the result and copy the device content result to it
    thrust::host_vector<CParticle3D_cuda> hostOutput;
    hostOutput = kernelOutput;

    // Print the results
    fprintf(stdOut, "\tPerform a modification (divide by 2, 4 and 8) of mass, m, position, x, and momentum p, of a list of atoms, using CUDA Thrust\r\n");
    for(size_t i=0; i<ARRAY_SIZE; i++)
    {
        fprintf(stdOut, "\t%i: m_in = %.2f; x_in = (%.2f, %.2f, %.2f); p_in = (%.2f, %.2f, %.2f); (m_in/2, x_in/4, p_in/8) = (%.2f, (%.2f, %.2f, %.2f), (%.2f, %.2f, %.2f))\r\n",
                int(i),
                hostAtomList[i].m_,
                hostAtomList[i].x_[0], hostAtomList[i].x_[1], hostAtomList[i].x_[2],
                hostAtomList[i].p_[0], hostAtomList[i].p_[1], hostAtomList[i].p_[2],
                hostOutput[i].m_,
                hostOutput[i].x_[0], hostOutput[i].x_[1], hostOutput[i].x_[2],
                hostOutput[i].p_[0], hostOutput[i].p_[1], hostOutput[i].p_[2]);
    }
}

int CCudaTest_cuda::getArraySize()
{
    return ARRAY_SIZE;
}
