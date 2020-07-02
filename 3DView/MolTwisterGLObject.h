#pragma once
#include <memory>
#include "Utilities/3DVector.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CGLObject
{
public:
    enum EType { objNone=0, objLine, objRect, objBox, objSphere };
    
public:
    CGLObject() { type_ = objNone; }

public:
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::shared_ptr<CGLObject> createCopy() const = 0;
    virtual std::shared_ptr<CGLObject> createCopy(const CGLObject& src) const = 0;

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
    virtual void serialize(CSerializer& io, bool saveToStream);
    virtual std::shared_ptr<CGLObject> createCopy() const { auto ret = std::shared_ptr<CGLObject>(new CGLObjectLine); *(CGLObjectLine*)(ret.get()) = *this; return ret; }
    virtual std::shared_ptr<CGLObject> createCopy(const CGLObject& src) const { auto ret = std::shared_ptr<CGLObject>(new CGLObjectLine(src)); *(CGLObjectLine*)(ret.get()) = *this; return ret; }

public:
    C3DVector p1_;
    C3DVector p2_;
    float color_[3];
};

END_CUDA_COMPATIBLE()
