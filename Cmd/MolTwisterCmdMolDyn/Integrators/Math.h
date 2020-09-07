#pragma once
#include <stdio.h>
#include <math.h>

#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMt
{
public:
    static double exp(double x);
    static double sinhXoverX(double x);
    
private:
    static double a_[6];
};

END_CUDA_COMPATIBLE()
