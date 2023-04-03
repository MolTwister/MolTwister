//
// Copyright (C) 2023 Richard Olsen.
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

#pragma once
#include <stdio.h>
#include <memory.h>
#ifdef __CUDACC__
    #include <cuda.h>
    #include <thrust/host_vector.h>
    #include <thrust/device_vector.h>
    #define mttransform thrust::transform
    #define mthost_vector thrust::host_vector
    #define mtdevice_vector thrust::device_vector

    #define mtraw_pointer_cast thrust::raw_pointer_cast

    #define mtcheck_env() printf("\tRunning as CU code...\r\n");

    #define mtblockDim blockDim
    #define mtblockIdx blockIdx
    #define mtthreadIdx threadIdx

    #define mtcudaMalloc cudaMalloc
    #define mtcudaMemcpy cudaMemcpy
    #define mtcudaFree cudaFree
    #define mtcudaDeviceSynchronize cudaDeviceSynchronize

    #define CUDA_PROC(n, m) <<<n, m>>>

    #define BEGIN_CUDA_COMPATIBLE() namespace mtdev {
    #define END_CUDA_COMPATIBLE() }

    #define HOSTDEV_CALLABLE __host__ __device__
    #define HOST_CALLABLE __host__
    #define DEV_CALLABLE __device__
    #define CUDA_GLOBAL __global__

    #define EXEC_POLICY thrust::device,
#else
    #include <algorithm>
    #include <vector>
    #define mttransform std::transform
    #define mthost_vector std::vector
    #define mtdevice_vector std::vector

    template<class T> T* mtraw_pointer_cast(const T* ptr) { return (T*)ptr; }

    #define mtcheck_env() printf("\tRunning as CPP code...\r\n");

    struct CCUDADim
    {
        CCUDADim() { x = y = z = 0; }
        long x, y, z;
    };
    #define mtblockDim CCUDADim()
    #define mtblockIdx CCUDADim()
    #define mtthreadIdx CCUDADim()
    #define cudaMemcpyDeviceToHost 0
    #define cudaMemcpyHostToHost 1
    #define cudaMemcpyHostToDevice 2
    #define cudaMemcpyDeviceToDevice 3

    #define mtcudaMalloc(ptr, size) { *(ptr) = malloc(size); }
    #define mtcudaMemcpy(dst, src, count, type) { memcpy(dst, src, count); }
    #define mtcudaFree free
    #define mtcudaDeviceSynchronize()

    #define CUDA_PROC(n, m)

    #define BEGIN_CUDA_COMPATIBLE()
    #define END_CUDA_COMPATIBLE()

    #define HOSTDEV_CALLABLE
    #define HOST_CALLABLE
    #define DEV_CALLABLE
    #define CUDA_GLOBAL

    #define EXEC_POLICY
#endif
