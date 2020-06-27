#include <iostream>
#include "MolTwisterGLObject.h"

BEGIN_CUDA_COMPATIBLE()

CGLObjectLine::CGLObjectLine(C3DVector p1, C3DVector p2, float r, float g, float b) : CGLObject()
{ 
    type_ = objLine;
    p1_ = p1; 
    p2_ = p2;
    
    color_[0] = r;
    color_[1] = g;
    color_[2] = b;
}

CGLObjectLine::CGLObjectLine(const CGLObject& src) : CGLObject()
{
    CGLObjectLine*  p;
    
    if(src.getType() == objLine)
    {
        p = (CGLObjectLine*)&src;
        
        type_ = objLine;
        *this = *p;
    }
}

END_CUDA_COMPATIBLE()
