#pragma once
#include <stdio.h>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class COut
{
public:
    static void printf(const char* format, ...);
    static void setOutputFile(FILE* outFile);
    static void setStdOut(FILE* stdOut);
    static void closeOutputFile() { if(outFile_) fclose(outFile_); }

private:
    static FILE* outFile_;
    static FILE* stdOut_;
};

END_CUDA_COMPATIBLE()
