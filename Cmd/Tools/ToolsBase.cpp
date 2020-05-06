#include "ToolsBase.h"

CToolsBase::CToolsBase(CMolTwisterState* state, FILE* stdOut)
{
    state_ = state;
    stdOut_ = stdOut;
}
