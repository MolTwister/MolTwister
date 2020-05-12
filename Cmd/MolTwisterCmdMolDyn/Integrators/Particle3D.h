#ifndef __ThesisMDTests__Particle3D__
#define __ThesisMDTests__Particle3D__

#include "../../../Utilities/3DVector.h"
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

class CParticle3D
{
public:
    CParticle3D();

public:
    double m;
    C3DVector x;
    C3DVector p;
};


#endif
