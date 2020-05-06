#include "3DBasis.h"

void C3DBasis::generateCartessianBasisAt(C3DVector, C3DVector newBasisZ)
{
    // Calculate Z unit vector of new basis, w, described in old system
    w_ = newBasisZ;
    w_.normalize();

    // Find a vector not parallel to NewBasisZ
    C3DVector l_p, l = newBasisZ, l_pxl;
    l_p = l + C3DVector(1.0, 0.0, 0.0);
    l_pxl = l_p.cross(l);
    if(l_pxl.isZero())
    {
        l_p = l + C3DVector(0.0, 1.0, 0.0);
        l_pxl = l_p.cross(l);
        if(l_pxl.isZero())
        {
            l_p = l + C3DVector(0.0, 0.0, 1.0);
            l_pxl = l_p.cross(l);
        }
    }

    // Calculate X unit vector of new basis, u, described in old system
    u_ = l_pxl;
    u_.normalize();

    // Calculate Y unit vector of new basis, v, described in old system
    v_ = w_.cross(u_);
}
