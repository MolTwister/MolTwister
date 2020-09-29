#include "Printf.h"
#include <stdio.h>
#include <cstdarg>

BEGIN_CUDA_COMPATIBLE()

FILE* COut::outFile_ = nullptr;
FILE* COut::stdOut_ = stdout;

void COut::printf(const char* format, ...)
{
    va_list     args;
    
    va_start(args, format);
    vfprintf(stdOut_, format, args);
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

void COut::setStdOut(FILE* stdOut)
{
    stdOut_ = stdOut;
}

END_CUDA_COMPATIBLE()
