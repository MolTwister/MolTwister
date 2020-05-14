#ifndef ThesisMDTests_Printf_h
#define ThesisMDTests_Printf_h

#include <stdio.h>

class COut
{
public:
    static void Printf(const char* format, ...);
    static void SetOutputFile(FILE* pOutFile);
    static void CloseOutputFile() { if(m_pOutFile) fclose(m_pOutFile); }

private:
    static FILE* m_pOutFile;
};

#endif
