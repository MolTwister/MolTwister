#pragma once
#include "MolTwisterState.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CToolsBase
{
public:
    CToolsBase() = delete;
    CToolsBase(CMolTwisterState* state, FILE* stdOut);

protected:
    CMolTwisterState* state_;
    FILE* stdOut_;
};

END_CUDA_COMPATIBLE()
