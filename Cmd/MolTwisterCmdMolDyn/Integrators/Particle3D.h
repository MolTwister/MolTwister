#pragma once
#include "../../../Utilities/3DVector.h"
#include "../../../Utilities/CUDAGeneralizations.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

BEGIN_CUDA_COMPATIBLE()

class CParticle3D
{
public:
    CParticle3D();

public:
    double m_;
    C3DVector x_;
    C3DVector p_;
};

END_CUDA_COMPATIBLE()
