#pragma once
#include <vector>
#include <string>
#include "CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CSerializer
{
public:
    CSerializer();

public:
    void operator<<(int i);
    void operator<<(size_t i);
    void operator<<(char c);
    void operator<<(double f);
    void operator<<(float f);
    void operator<<(bool b);
    void operator<<(std::string s);

    void operator>>(int& i);
    void operator>>(size_t& i);
    void operator>>(char& c);
    void operator>>(double& f);
    void operator>>(float& f);
    void operator>>(bool& b);
    void operator>>(std::string& s);

    void seekBegin();

private:
    std::vector<std::string> archive_;
    int readPos_;
};

END_CUDA_COMPATIBLE()
