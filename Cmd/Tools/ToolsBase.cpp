#include "ToolsBase.h"

BEGIN_CUDA_COMPATIBLE()

CToolsBase::CToolsBase(CMolTwisterState* state, FILE* stdOut)
{
    state_ = state;
    stdOut_ = stdOut;
}

END_CUDA_COMPATIBLE()
