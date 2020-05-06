#pragma once
#include <vector>

class CExpLookup
{
public:
    CExpLookup() { delta_ = 1.0; low_ = 0.0; }
    
public:
    void init(double lowArg, double highArg, int N);
    double exp(double arg) const;
    
private:
    double delta_;
    double low_;
    std::vector<double> values_;
};
