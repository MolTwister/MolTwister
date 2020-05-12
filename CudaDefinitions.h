#pragma once

#ifdef __CUDACC__
#include <cuda.h>
#define HOSTDEV_CALLABLE __host__ __device__
#define HOST_CALLABLE __host__
#define DEV_CALLABLE __device__
#define CUDA_GLOBAL __global__
#else
#define HOSTDEV_CALLABLE
#define HOST_CALLABLE
#define DEV_CALLABLE
#define CUDA_GLOBAL
#endif
