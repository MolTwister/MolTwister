#include <iostream>
#include <math.h>
#include "ExpLookup.h"

void CExpLookup::init(double lowArg, double highArg, int N)
{
    if(N == 0) return;
    
    values_.reserve(N);
    
    delta_ = (highArg - lowArg) / double(N);
    low_ = lowArg;
    
    for(int i=0; i<N; i++)
    {
        values_.emplace_back(::exp(lowArg + delta_ * double(i)));
    }
}

double CExpLookup::exp(double arg) const
{
    int n = int((arg - low_) / delta_);

    if((n < 0) || (n >= values_.size())) return ::exp(arg);
    
    return values_[n];
}
