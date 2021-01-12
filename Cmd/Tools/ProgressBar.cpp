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

#include "ProgressBar.h"

void CProgressBar::beginProgress(std::string descript)
{
    printf("\r\n\t %s:\r\n", descript.data());
    lastNumDots_ = -1;
}

void CProgressBar::updateProgress(int step, int totSteps)
{
    const int maxDots = 50;
    double fracDone = double(step) / double(totSteps);
    int numDots = (int)(double(maxDots) * fracDone);

    if(numDots != lastNumDots_)
    {
        printf("\r\t|");
        for(int i=0; i<maxDots; i++)
        {
            if(i < numDots) printf("*");
            else             printf(" ");
        }
        printf("| [%.3i%%]", (int)(fracDone*100.0));
        fflush(stdout);

        lastNumDots_ = numDots;
    }
}

void CProgressBar::endProgress()
{
    printf("\r\t|**************************************************| [100%%]\r\n");
}
