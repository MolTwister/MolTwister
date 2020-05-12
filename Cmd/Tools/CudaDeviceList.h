#pragma once
#include <stdio.h>
#include <vector>
#include <string>
#include <memory>

class CCudaDeviceInfo
{
public:
    CCudaDeviceInfo();
    CCudaDeviceInfo(int iDevice);

public:
    void Print(FILE* pFile=stdout) const;
    std::string	GetName() const { return m_szName; }
    int	GetCudaDevice() const { return m_iCudaDevice; }
    int	GetCudaRuntimeVersion() const { return m_iCudaRuntimeVersion; }
    int	GetStackSize() const { return m_iStackSize; }
    int	GetPrintfFIFOSize() const { return m_iPrintfFIFOSize; }
    int	GetMallocHeapSize() const { return m_iMallocHeapSize; }
    int	GetTotalMem() const { return m_iTotalMem; }
    int	GetFreeMem() const { return m_iFreeMem; }
    int	GetMaxSharedMemPerBlock() const { return m_iMaxSharedMemPerBlock; }
    int	GetMaxRegsPerBlock() const { return m_iMaxRegsPerBlock; }
    int	GetThreadsPerBlock() const { return m_iThreadsPerBlock; }
    int	GetMaxThreadsDimX() const { return m_iMaxThreadsDim[0]; }
    int	GetMaxThreadsDimY() const { return m_iMaxThreadsDim[1]; }
    int	GetMaxThreadsDimZ() const { return m_iMaxThreadsDim[2]; }
    int	GetMaxGridDimX() const { return m_iMaxGridDim[0]; }
    int	GetMaxGridDimY() const { return m_iMaxGridDim[1]; }
    int	GetMaxGridDimZ() const { return m_iMaxGridDim[2]; }
    int	GetPeakClockFreq_kHz() const { return m_iPeakClockFreq_kHz; }
    int	GetMaxConstantMem() const { return m_iMaxConstantMem; }
    int	GetNumMultiProcsOnDevice() const { return m_iNumMultiProcsOnDevice; }
    int	GetPCIBusID() const { return m_iPCIBusID; }
    int	GetPCIDevID() const { return m_iPCIDevID; }
    int	GetPeakMemClockFreq_kHz() const { return m_iPeakMemClockFreq_kHz; }
    int	GetMemBusWidth() const { return m_iMemBusWidth; }
    int	GetMaxThreadsPerMultiProc() const { return m_iMaxThreadsPerMultiProc; }

private:
    void ResetToDefault();

private:
    std::string	m_szName;
    int	m_iCudaDevice;
    int	m_iCudaRuntimeVersion;
    int	m_iStackSize;
    int	m_iPrintfFIFOSize;
    int	m_iMallocHeapSize;
    int	m_iTotalMem;
    int	m_iFreeMem;
    int	m_iMaxSharedMemPerBlock;
    int	m_iMaxRegsPerBlock;
    int	m_iThreadsPerBlock;
    int	m_iMaxThreadsDim[3];
    int	m_iMaxGridDim[3];
    int	m_iPeakClockFreq_kHz;
    int	m_iMaxConstantMem;
    int	m_iNumMultiProcsOnDevice;
    int	m_iPCIBusID;
    int	m_iPCIDevID;
    int	m_iPeakMemClockFreq_kHz;
    int	m_iMemBusWidth;
    int	m_iMaxThreadsPerMultiProc;
};


class CCudaDeviceList
{
public:
    CCudaDeviceList();

public:
    int GetActiveDevice() const;
    static bool SetActiveDevice(int iIndex);
    int GetNumDevices() const;
    CCudaDeviceInfo& GetDevice(int iIndex) const;
    CCudaDeviceInfo& operator[](int iIndex) const;
    void ResetActiveGPU();
    
    long GetThreadsPerBlock();
    long GetRegistersPerBlock();
    std::pair<long, long> CalcNumGPUBlocks1D(long lNumItems, short nRegisterSizeLimit=-1);

private:
    std::vector<std::shared_ptr<CCudaDeviceInfo>> m_aDevices;
};

