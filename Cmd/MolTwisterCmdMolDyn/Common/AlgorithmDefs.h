#pragma once
#include <stdio.h>
#ifdef __CUDACC__
    #include <cuda.h>
    #include <thrust/host_vector.h>
    #include <thrust/device_vector.h>
    #define mttransform thrust::transform
    #define mthost_vector thrust::host_vector
    #define mtdevice_vector thrust::device_vector

    #define mtraw_pointer_cast thrust::raw_pointer_cast

    #define mtcheck_env() printf("\tRunning as CU code...\r\n");

    #define HOSTDEV_CALLABLE __host__ __device__
    #define HOST_CALLABLE __host__
    #define DEV_CALLABLE __device__
    #define CUDA_GLOBAL __global__
#else
    #include <algorithm>
    #include <vector>
    #define mttransform std::transform
    #define mthost_vector std::vector
    #define mtdevice_vector std::vector

    template<class T> T* mtraw_pointer_cast(const T* ptr) { return (T*)ptr; }

    #define mtcheck_env() printf("\tRunning as CPP code...\r\n");

    #define HOSTDEV_CALLABLE
    #define HOST_CALLABLE
    #define DEV_CALLABLE
    #define CUDA_GLOBAL
#endif
