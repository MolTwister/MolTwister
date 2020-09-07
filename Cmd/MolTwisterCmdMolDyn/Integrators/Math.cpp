#include "Math.h"

#define EXP1 2.71828

BEGIN_CUDA_COMPATIBLE()

double CMt::a_[6] = { 1.0, 1.0/6.0, 1.0/120.0, 1.0/5040.0, 1.0/362880.0, 1.0/39916800.0 };

double CMt::exp(double x)
{
    // Improve numerical stability by linear approximation
    // of exp(x), around x=1, for x>1. Valid for equations
    // where x is small during natural progression.
    if(x > 1.0)
        return EXP1 * x;
    
    return ::exp(x);
}

double CMt::sinhXoverX(double x)
{
    double sum = 0.0;
    double x_n = 1.0;
    
    // Improve numerical stability by rather calculating
    // (1/x)*TaylorExpansion(sinh(x)) to 10th order
    for(int n=0; n<6; n++)
    {
        sum+= a_[n]*x_n;
        x_n*= (x*x);
    }
    
    return sum;
}

END_CUDA_COMPATIBLE()
