// Note! This is an entry point to GPU compiled code.
// Hence, it is necessary to include this file twice
// in different namespaces, within the same file. This
// means that 'pragma once' should not be used

#include <string>
#include "../../../Utilities/CUDAGeneralizations.h"
#include "../../../Utilities/Serializer.h"
#include "../Config/MolDynConfigStruct.h"

BEGIN_CUDA_COMPATIBLE()

class CMDSimulator
{
public:
    CMDSimulator() {}

public:
    static void run(SMolDynConfigStruct config, FILE* stdOut, CSerializer& stateContent);
    static void run(SMolDynConfigStruct config, FILE* stdOut, void* state);
};

END_CUDA_COMPATIBLE()
