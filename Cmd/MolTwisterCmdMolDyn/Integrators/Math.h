#pragma once
#include <stdio.h>
#include <math.h>

#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMt
{
public:
    static double Exp(double x);
    static double SinhXoverX(double x);
    
private:
    static double a[6];
};

END_CUDA_COMPATIBLE()
