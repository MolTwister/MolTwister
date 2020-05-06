#pragma once
#include "MolTwisterState.h"

class CToolsBase
{
public:
    CToolsBase() = delete;
    CToolsBase(CMolTwisterState* state, FILE* stdOut);

protected:
    CMolTwisterState* state_;
    FILE* stdOut_;
};
