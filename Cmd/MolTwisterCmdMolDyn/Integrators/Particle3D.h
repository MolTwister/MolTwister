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
    double m;
    C3DVector x;
    C3DVector p;
};

END_CUDA_COMPATIBLE()
