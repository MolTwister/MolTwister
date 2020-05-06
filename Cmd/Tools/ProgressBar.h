#pragma once
#include <string>

class CProgressBar
{
public:
    CProgressBar() { lastNumDots_ = -1; }

public:
    void beginProgress(std::string descript);
    void updateProgress(int step, int totSteps);
    void endProgress();

private:
    int lastNumDots_;
};
