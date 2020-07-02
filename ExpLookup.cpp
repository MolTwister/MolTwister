#include <iostream>
#include <math.h>
#include "ExpLookup.h"

void CExpLookup::serialize(std::stringstream& io, bool saveToStream)
{
    if(saveToStream)
    {
        io << values_.size();
        for(double val : values_)
        {
            io << val;
        }

        io << delta_;
        io << low_;
    }
    else
    {
        size_t size;
        io >> size;
        values_.resize(size);
        for(size_t i=0; i<size; i++)
        {
            double val;
            io >> val;
            values_[i] = val;
        }

        io >> delta_;
        io >> low_;
    }
}

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

    if((n < 0) || (n >= (int)values_.size())) return ::exp(arg);
    
    return values_[n];
}
