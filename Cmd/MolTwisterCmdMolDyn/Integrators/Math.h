#pragma once
#include <stdio.h>
#include <math.h>

class CMt
{
public:
    static double Exp(double x);
    static double SinhXoverX(double x);
    
private:
    static double a[6];
};
