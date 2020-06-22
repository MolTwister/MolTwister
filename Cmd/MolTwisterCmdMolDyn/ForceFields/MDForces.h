#pragma once
#include "../../../Utilities/3DVector.h"

class CMDForces
{
public:
    CMDForces() {}
    CMDForces(const C3DVector& F, const C3DVector& Fpi) { F_ = F; Fpi_ = Fpi; }

public:
    C3DVector F_;
    C3DVector Fpi_;
};
