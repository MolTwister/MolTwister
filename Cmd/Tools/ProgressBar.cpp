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
