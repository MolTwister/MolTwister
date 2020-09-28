#include "Printf.h"
#include <stdio.h>
#include <cstdarg>

BEGIN_CUDA_COMPATIBLE()

FILE* COut::outFile_ = nullptr;

void COut::printf(const char* format, ...)
{
    va_list     args;
    
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    if(outFile_)
    {
        va_start(args, format);
        vfprintf(outFile_, format, args);
        va_end(args);
    }
    fflush(outFile_);
}

void COut::setOutputFile(FILE* outFile)
{
    if(!outFile) printf("Warning! Could not open output file!\r\n");
    outFile_ = outFile;
}

END_CUDA_COMPATIBLE()
