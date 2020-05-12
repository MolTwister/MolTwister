#include <cuda.h>
#include <utility>
#include "CudaDeviceList.h"

CCudaDeviceInfo::CCudaDeviceInfo()
{
    ResetToDefault();
}

void CCudaDeviceInfo::ResetToDefault()
{
    m_szName = "";
    m_iCudaDevice = -1;
    m_iCudaRuntimeVersion = -1;
    m_iStackSize = -1;
    m_iPrintfFIFOSize = -1;
    m_iMallocHeapSize = -1;
    m_iTotalMem = -1;
    m_iFreeMem = -1;
    m_iMaxSharedMemPerBlock = -1;
    m_iMaxRegsPerBlock = -1;
    m_iThreadsPerBlock = -1;
    for(int i=0; i<3; i++) m_iMaxThreadsDim[i] = -1;
    for(int i=0; i<3; i++) m_iMaxGridDim[i] = -1;
    m_iPeakClockFreq_kHz = -1;
    m_iMaxConstantMem = -1;
    m_iNumMultiProcsOnDevice = -1;
    m_iPCIBusID = -1;
    m_iPCIDevID = -1;
    m_iPeakMemClockFreq_kHz = -1;
    m_iMemBusWidth = -1;
    m_iMaxThreadsPerMultiProc = -1;
}

CCudaDeviceInfo::CCudaDeviceInfo(int iDevice)
{
    cudaDeviceProp  DevProp;
    cudaError       ErrorCode;
    size_t          Size1, Size2;


    ResetToDefault();

    if(iDevice == -1)
    {
        cudaGetDevice(&iDevice);
    }

    m_iCudaDevice = iDevice;
    cudaRuntimeGetVersion(&m_iCudaRuntimeVersion);

    ErrorCode = cudaDeviceGetLimit(&Size1, cudaLimitStackSize);
    if(ErrorCode == cudaSuccess) m_iStackSize = (int)Size1;

    ErrorCode = cudaDeviceGetLimit(&Size1, cudaLimitPrintfFifoSize);
    if(ErrorCode == cudaSuccess) m_iPrintfFIFOSize = (int)Size1;

    ErrorCode = cudaDeviceGetLimit(&Size1, cudaLimitMallocHeapSize);
    if(ErrorCode == cudaSuccess) m_iMallocHeapSize = (int)Size1;

    ErrorCode = cudaMemGetInfo(&Size1, &Size2);
    if(ErrorCode == cudaSuccess) { m_iFreeMem = (int)Size1; m_iTotalMem = (int)Size2; }

    ErrorCode = cudaGetDeviceProperties(&DevProp, iDevice);
    if(ErrorCode == cudaSuccess)
    {
        char szName[257];
        for(int i=0; i<256; i++) szName[i] = DevProp.name[i];
        szName[256] = '\0';
        m_szName = szName;
        m_iMaxSharedMemPerBlock = (int)DevProp.sharedMemPerBlock;
        m_iMaxRegsPerBlock = DevProp.regsPerBlock;
        m_iThreadsPerBlock = DevProp.maxThreadsPerBlock;
        for(int i=0; i<3; i++) m_iMaxThreadsDim[i] = DevProp.maxThreadsDim[i];
        for(int i=0; i<3; i++) m_iMaxGridDim[i] = DevProp.maxGridSize[i];
        m_iPeakClockFreq_kHz = DevProp.clockRate;
        m_iMaxConstantMem = DevProp.totalConstMem;
        m_iNumMultiProcsOnDevice = DevProp.multiProcessorCount;
        m_iPCIBusID = DevProp.pciBusID;
        m_iPCIDevID = DevProp.pciDeviceID;
        m_iPeakMemClockFreq_kHz = DevProp.memoryClockRate;
        m_iMemBusWidth = DevProp.memoryBusWidth;
        m_iMaxThreadsPerMultiProc = DevProp.maxThreadsPerMultiProcessor;
    }
}

void CCudaDeviceInfo::Print(FILE* pFile) const
{
    fprintf(pFile, "\tCUDA device limits: Ver %.4f, Device: %i\r\n", double(m_iCudaRuntimeVersion) / 1000.0, m_iCudaDevice);
    fprintf(pFile, "\t---------------------------------------------\r\n");

    fprintf(pFile, "\tDevice name: %s\r\n", m_szName.data());

    if(m_iPCIBusID != -1)
        fprintf(pFile, "\tPCI bus ID: %i\r\n", m_iPCIBusID);

    if(m_iPCIDevID != -1)
        fprintf(pFile, "\tPCI device ID: %i\r\n", m_iPCIDevID);


    if(m_iStackSize != -1)
        fprintf(pFile, "\tStack size: %i bytes, %.2f kb, %.2f MB\r\n", (int)m_iStackSize, double(m_iStackSize)/1024.0, double(m_iStackSize)/(1024.0*1024.0));

    if(m_iPrintfFIFOSize != -1)
        fprintf(pFile, "\tPrintf-FIFO size: %i, %.2f kb, %.2f MB\r\n", (int)m_iPrintfFIFOSize, double(m_iPrintfFIFOSize)/1024.0, double(m_iPrintfFIFOSize)/(1024.0*1024.0));

    if(m_iMallocHeapSize != -1)
        fprintf(pFile, "\tMalloc heap size: %i, %.2f kb, %.2f MB\r\n", (int)m_iMallocHeapSize, double(m_iMallocHeapSize)/1024.0, double(m_iMallocHeapSize)/(1024.0*1024.0));

    if(m_iTotalMem != -1)
        fprintf(pFile, "\tTotal memory: %i, %.2f kb, %.2f MB\r\n", (int)m_iTotalMem, double(m_iTotalMem)/1024.0, double(m_iTotalMem)/(1024.0*1024.0));

    if(m_iFreeMem != -1)
        fprintf(pFile, "\tFree memory: %i, %.2f kb, %.2f MB\r\n", (int)m_iFreeMem, double(m_iFreeMem)/1024.0, double(m_iFreeMem)/(1024.0*1024.0));

    if(m_iMaxConstantMem != -1)
        fprintf(pFile, "\tConstant memory: %i, %.2f kb, %.2f MB\r\n", (int)m_iMaxConstantMem, double(m_iMaxConstantMem)/1024.0, double(m_iMaxConstantMem)/(1024.0*1024.0));

    if(m_iMemBusWidth != -1)
        fprintf(pFile, "\tMemory bus width: %i bits\r\n", m_iMemBusWidth);


    if(m_iMaxSharedMemPerBlock != -1)
        fprintf(pFile, "\tMax shared memory per block: %i, %.2f kb, %.2f MB\r\n", (int)m_iMaxSharedMemPerBlock, double(m_iMaxSharedMemPerBlock)/1024.0, double(m_iMaxSharedMemPerBlock)/(1024.0*1024.0));

    if(m_iMaxRegsPerBlock != -1)
        fprintf(pFile, "\tMaximium registers per block: %i\r\n", m_iMaxRegsPerBlock);

    if(m_iThreadsPerBlock != -1)
        fprintf(pFile, "\tThreads per block: %i\r\n", m_iThreadsPerBlock);

    if(m_iMaxThreadsDim[0] != -1)
        fprintf(pFile, "\tThread dimensions: (%i, %i, %i)\r\n", m_iMaxThreadsDim[0], m_iMaxThreadsDim[1], m_iMaxThreadsDim[2]);

    if(m_iMaxGridDim[0] != -1)
        fprintf(pFile, "\tGrid dimensions: (%i, %i, %i)\r\n", m_iMaxGridDim[0], m_iMaxGridDim[1], m_iMaxGridDim[2]);


    if(m_iPeakClockFreq_kHz != -1)
        fprintf(pFile, "\tPeak clock frequency: %i kHz\r\n", m_iPeakClockFreq_kHz);

    if(m_iPeakMemClockFreq_kHz != -1)
        fprintf(pFile, "\tPeak memory clock frequency: %i kHz\r\n", m_iPeakMemClockFreq_kHz);


    if(m_iNumMultiProcsOnDevice != -1)
        fprintf(pFile, "\tNum. of multiprocs. on device: %i\r\n", m_iNumMultiProcsOnDevice);

    if(m_iMaxThreadsPerMultiProc != -1)
        fprintf(pFile, "\tMax threads on one multiproc.: %i\r\n", m_iMaxThreadsPerMultiProc);
}



CCudaDeviceList::CCudaDeviceList()
{
    int iNumDevices;

    cudaGetDeviceCount(&iNumDevices);
    
    for(int i=0; i<iNumDevices; i++)
    {
        std::shared_ptr<CCudaDeviceInfo> CudaDeviceInfo = std::shared_ptr<CCudaDeviceInfo>(new CCudaDeviceInfo(i));

        m_aDevices.push_back(CudaDeviceInfo);
    }
}

int CCudaDeviceList::GetActiveDevice() const
{
    int iDevice = -1;

    cudaGetDevice(&iDevice);

    return iDevice;
}

bool CCudaDeviceList::SetActiveDevice(int iIndex)
{
    if(cudaSetDevice(iIndex) != cudaSuccess) return false;
    return true;
}

int CCudaDeviceList::GetNumDevices() const
{
    return static_cast<int>(m_aDevices.size());
}

CCudaDeviceInfo& CCudaDeviceList::GetDevice(int iIndex) const
{
    return *m_aDevices[iIndex];
}

CCudaDeviceInfo& CCudaDeviceList::operator[](int iIndex) const
{
    return *m_aDevices[iIndex];
}

void CCudaDeviceList::ResetActiveGPU()
{
    int iDevice = -1;

    cudaGetDevice(&iDevice);
    if(cudaDeviceReset() == cudaSuccess)
        printf("\tSuccesfully reset GPU %i\r\n", iDevice);
    else
        printf("\tFailed to reset GPU %i\r\n", iDevice);
}

long CCudaDeviceList::GetThreadsPerBlock()
{
    int iActiveDev = GetActiveDevice();
    if((iActiveDev < 0) || (iActiveDev >= GetNumDevices()))
    {
        printf("Could not find an active CUDA device");
        return -1;
    }

    long lThreadsPerBlock = static_cast<long>(GetDevice(iActiveDev).GetThreadsPerBlock());
    if(lThreadsPerBlock <= 0)
    {
        printf("Found no CUDA threads per CUDA block");
        return -1;
    }

    return lThreadsPerBlock;
}

long CCudaDeviceList::GetRegistersPerBlock()
{
    int iActiveDev = GetActiveDevice();
    if((iActiveDev < 0) || (iActiveDev >= GetNumDevices()))
    {
        printf("Could not find an active CUDA device");
        return -1;
    }

    long lRegistersPerBlock = static_cast<long>(GetDevice(iActiveDev).GetMaxRegsPerBlock());
    if(lRegistersPerBlock <= 0)
    {
        printf("Found no CUDA registers per CUDA block");
        return -1;
    }

    return lRegistersPerBlock;
}

std::pair<long, long> CCudaDeviceList::CalcNumGPUBlocks1D(long lNumItems, short nRegisterSizeLimit)
{
    // Retrieve the number of threads per block
    long lThreadsPerBlock;
    if(nRegisterSizeLimit < 0)
        lThreadsPerBlock = GetThreadsPerBlock();
    else
    {
        long lThreadsPerBlockHW = GetThreadsPerBlock();
        lThreadsPerBlock = (long)ceil(double(GetRegistersPerBlock()) / double(nRegisterSizeLimit));
        if(lThreadsPerBlock > lThreadsPerBlockHW) lThreadsPerBlock = lThreadsPerBlockHW;
    }

    // Calculate the number of blocks. Should always overestimate!
    long lNumBlocks = lNumItems / lThreadsPerBlock + ((lNumItems % lThreadsPerBlock) ? 1 : 0);

    return std::pair<long, long>(lNumBlocks, lThreadsPerBlock);
}
