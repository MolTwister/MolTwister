#pragma once
#include <vector>
#include "Utilities/Serializer.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CExpLookup
{
public:
    CExpLookup() { delta_ = 1.0; low_ = 0.0; }
    
public:
    void serialize(CSerializer& io, bool saveToStream);
    void init(double lowArg, double highArg, int N);
    double exp(double arg) const;
    
private:
    double delta_;
    double low_;
    std::vector<double> values_;
};

END_CUDA_COMPATIBLE()
