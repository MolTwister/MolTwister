//
// Copyright (C) 2021 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

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
