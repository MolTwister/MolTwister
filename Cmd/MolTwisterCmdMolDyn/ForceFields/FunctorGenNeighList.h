#pragma once
#include "../../../Utilities/CUDAGeneralizations.h"
#include "MDFFMatrices.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenNeighList
{
public:
    HOSTDEV_CALLABLE CFunctorGenNeighList() {}

public:
    HOSTDEV_CALLABLE int operator()(CMDFFMatrices::CAtom& atom);
};

END_CUDA_COMPATIBLE()
