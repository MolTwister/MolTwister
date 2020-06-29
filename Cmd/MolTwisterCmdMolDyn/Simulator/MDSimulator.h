// Note! This is an entry point to GPU compiled code.
// Hence, it is necessary to include this file twice
// in different namespaces, within the same file. This
// means that 'pragma once' should not be used

#include <string>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMDSimulator
{
public:
    enum Ensamble { ensambleNVT=0, ensambleNPT=1 };

public:
    CMDSimulator() {}

public:
    static void run(int numMDSteps, int outputStride, int nhChainLength, std::string progressOutFileName, FILE* stdOut, void* state, int ensamble);
};

END_CUDA_COMPATIBLE()
