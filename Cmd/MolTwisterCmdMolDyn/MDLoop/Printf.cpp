#include "Printf.h"
#include <stdio.h>
#include <cstdarg>

BEGIN_CUDA_COMPATIBLE()

FILE* COut::m_pOutFile = NULL;

void COut::Printf(const char* format, ...)
{
    va_list     args;
    
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    if(m_pOutFile)
    {
        va_start(args, format);
        vfprintf(m_pOutFile, format, args);
        va_end(args);
    }
    fflush(m_pOutFile);
}

void COut::SetOutputFile(FILE* pOutFile)
{
    if(!pOutFile) printf("Warning! Could not open output file!\r\n");
    m_pOutFile = pOutFile;
}

END_CUDA_COMPATIBLE()
