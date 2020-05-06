#pragma once
#include "Utilities/3DVector.h"

class CGLObject
{
public:
    enum EType { objNone=0, objLine, objRect, objBox, objSphere };
    
public:
    CGLObject() { type_ = objNone; }
    
public:
    EType getType() const { return type_; }
    
protected:
    EType type_;
};

class CGLObjectLine : public CGLObject
{
public:
    CGLObjectLine() : CGLObject() { type_ = objLine; }
    CGLObjectLine(C3DVector p1, C3DVector p2, float r, float g, float b);
    CGLObjectLine(const CGLObject& src);
    
public:
    C3DVector p1_;
    C3DVector p2_;
    float color_[3];
};
