#pragma once
#include "../../../Utilities/3DVector.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMDForces
{
public:
    HOSTDEV_CALLABLE CMDForces() {}
    HOSTDEV_CALLABLE CMDForces(const C3DVector& F, const C3DVector& Fpi) { F_ = F; Fpi_ = Fpi; }

public:
    C3DVector F_;
    C3DVector Fpi_;
};

END_CUDA_COMPATIBLE()
