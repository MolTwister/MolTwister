#pragma once
#include <stdio.h>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class COut
{
public:
    static void Printf(const char* format, ...);
    static void SetOutputFile(FILE* pOutFile);
    static void CloseOutputFile() { if(m_pOutFile) fclose(m_pOutFile); }

private:
    static FILE* m_pOutFile;
};

END_CUDA_COMPATIBLE()
