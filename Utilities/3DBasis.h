#pragma once
#include "3DVector.h"

class C3DBasis
{
public:
    C3DBasis() {}
    C3DBasis(C3DVector u, C3DVector v, C3DVector w) { u_ = u; v_ = v; w_ = w; }

public:
    void generateCartessianBasisAt(C3DVector pos, C3DVector newBasisZ);

public:
    C3DVector u_;
    C3DVector v_;
    C3DVector w_;
};
