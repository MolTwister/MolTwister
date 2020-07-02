#include <iostream>
#include <stdio.h>
#include <math.h>
#include "3DVector.h"
#include "3DRect.h"

BEGIN_CUDA_COMPATIBLE()

void C3DVector::print() const
{
    printf("(%.g, %.g, %.g)", x_, y_, z_);
}

HOSTDEV_CALLABLE double C3DVector::norm() const
{
    return sqrt(norm2());
}

HOSTDEV_CALLABLE double C3DVector::norm2() const
{
    return x_*x_ + y_*y_ + z_*z_;
}

HOSTDEV_CALLABLE double C3DVector::distToAcrossPBC(const C3DVector& v, const C3DRect& pbc) const
{
    double dX = (x_ - v.x_);
    double dY = (y_ - v.y_);
    double dZ = (z_ - v.z_);
    double w;
    
    
    dX = (dX < 0.0) ? -dX : dX;
    dY = (dY < 0.0) ? -dY : dY;
    dZ = (dZ < 0.0) ? -dZ : dZ;
    
    w = pbc.getWidthX();
    if(dX > (w / 2.0)) dX = w - dX;
    
    w = pbc.getWidthY();
    if(dY > (w / 2.0)) dY = w - dY;
    
    w = pbc.getWidthZ();
    if(dZ > (w / 2.0)) dZ = w - dZ;
    
    return sqrt(dX*dX + dY*dY + dZ*dZ);
}

HOSTDEV_CALLABLE double C3DVector::distToAcrossPBC2(const C3DVector& v, const C3DRect& pbc) const
{
    double dX = (x_ - v.x_);
    double dY = (y_ - v.y_);
    double dZ = (z_ - v.z_);
    double w;
    
    
    dX = (dX < 0.0) ? -dX : dX;
    dY = (dY < 0.0) ? -dY : dY;
    dZ = (dZ < 0.0) ? -dZ : dZ;
    
    w = pbc.getWidthX();
    if(dX > (w / 2.0)) dX = w - dX;
    
    w = pbc.getWidthY();
    if(dY > (w / 2.0)) dY = w - dY;
    
    w = pbc.getWidthZ();
    if(dZ > (w / 2.0)) dZ = w - dZ;
    
    return dX*dX + dY*dY + dZ*dZ;
}

HOSTDEV_CALLABLE double C3DVector::distTo(const C3DVector& v) const
{
    double dX = (x_ - v.x_);
    double dY = (y_ - v.y_);
    double dZ = (z_ - v.z_);

    return sqrt(dX*dX + dY*dY + dZ*dZ);
}

HOSTDEV_CALLABLE double C3DVector::cosAngle() const
{
    C3DVector v2(1.0, 0.0, 0.0);
    
    return cosAngle(v2);
}

HOSTDEV_CALLABLE double C3DVector::cosAngle(const C3DVector& v2) const
{
    double d = norm() * v2.norm();
    
    if(d == 0.0) return static_cast<double>(NAN);

    C3DVector thisVec(x_, y_, z_);
    return (thisVec * v2) / d;
}

HOSTDEV_CALLABLE double C3DVector::angle() const
{
    return acos(cosAngle());
}

HOSTDEV_CALLABLE double C3DVector::angle(const C3DVector& v2) const
{
    double cosAlpha = cosAngle(v2);
    
    if(cosAlpha > 1.0) cosAlpha = 1.0;
    if(cosAlpha < -1.0) cosAlpha = -1.0;
    
    return acos(cosAlpha);
}

HOSTDEV_CALLABLE void C3DVector::moveToSameSideOfPBCAsThis(C3DVector& v, const C3DRect& pbc) const
{
    double dX = (x_ - v.x_);
    double dY = (y_ - v.y_);
    double dZ = (z_ - v.z_);
    double w;
    
    
    dX = (dX < 0.0) ? -dX : dX;
    dY = (dY < 0.0) ? -dY : dY;
    dZ = (dZ < 0.0) ? -dZ : dZ;
    
    w = pbc.getWidthX();
    if(dX > (w / 2.0))
    {
        if(x_ > v.x_) v.x_+= w;
        else        v.x_-= w;
    }
    
    w = pbc.getWidthY();
    if(dY > (w / 2.0))
    {
        if(y_ > v.y_) v.y_+= w;
        else        v.y_-= w;
    }
    
    w = pbc.getWidthZ();
    if(dZ > (w / 2.0))
    {
        if(z_ > v.z_) v.z_+= w;
        else        v.z_-= w;
    }
}

HOSTDEV_CALLABLE double C3DVector::cosAngleAcrossPBC(C3DVector v1, C3DVector v3, const C3DRect &pbc) const
{
    moveToSameSideOfPBCAsThis(v1, pbc);
    moveToSameSideOfPBCAsThis(v3, pbc);
    
    C3DVector v2 = *this;
    C3DVector v12 = v1 - v2;
    C3DVector v32 = v3 - v2;
    
    return (v12 * v32) / ( v12.norm() * v32.norm() );
}

HOSTDEV_CALLABLE double C3DVector::angleAcrossPBC(const C3DVector &v1, const C3DVector &v3, const C3DRect &pbc) const
{
    double cosAlpha = cosAngleAcrossPBC(v1, v3, pbc);
    
    if(cosAlpha > 1.0) cosAlpha = 1.0;
    if(cosAlpha < -1.0) cosAlpha = -1.0;
    
    return acos(cosAlpha);
}

HOSTDEV_CALLABLE C3DVector C3DVector::unit() const
{
    const double nan = static_cast<double>(NAN);
    C3DVector ret(nan, nan, nan);

    double d = norm();
    if(d == 0.0) return ret;

    ret.x_ = x_ / d;
    ret.y_ = y_ / d;
    ret.z_ = z_ / d;

    return ret;
}

HOSTDEV_CALLABLE void C3DVector::normalize()
{
    double d = norm();
    if(d == 0.0) return;
    
    x_/= d;
    y_/= d;
    z_/= d;
}

HOSTDEV_CALLABLE C3DVector C3DVector::cross(const C3DVector &rhs) const
{
    C3DVector ret;
    
    ret.x_ = (y_*rhs.z_ - z_*rhs.y_);
    ret.y_ = (z_*rhs.x_ - x_*rhs.z_);
    ret.z_ = (x_*rhs.y_ - y_*rhs.x_);
    
    return ret;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator+(const C3DVector& rhs) const
{
    C3DVector ret;
    
    ret.x_ = x_ + rhs.x_;
    ret.y_ = y_ + rhs.y_;
    ret.z_ = z_ + rhs.z_;
    
    return ret;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator+=(const C3DVector& rhs)
{
    C3DVector ret;
    
    x_+= rhs.x_;
    y_+= rhs.y_;
    z_+= rhs.z_;
    
    return *this;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator-(const C3DVector& rhs) const
{
    C3DVector ret;
    
    ret.x_ = x_ - rhs.x_;
    ret.y_ = y_ - rhs.y_;
    ret.z_ = z_ - rhs.z_;
    
    return ret;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator-=(const C3DVector& rhs)
{
    C3DVector ret;
    
    x_-= rhs.x_;
    y_-= rhs.y_;
    z_-= rhs.z_;
    
    return *this;
}

HOSTDEV_CALLABLE double C3DVector::operator*(const C3DVector& rhs) const
{
    return (x_*rhs.x_) + (y_*rhs.y_) + (z_*rhs.z_);
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator*(double a) const
{
    C3DVector ret;
    
    ret.x_ = x_ * a;
    ret.y_ = y_ * a;
    ret.z_ = z_ * a;
    
    return ret;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator*=(double a)
{
    x_*= a;
    y_*= a;
    z_*= a;
    
    return *this;
}

HOSTDEV_CALLABLE bool C3DVector::operator==(const C3DVector& rhs) const
{
    if(fabs(x_ - rhs.x_) > cmpLimit_) return false;
    if(fabs(y_ - rhs.y_) > cmpLimit_) return false;
    if(fabs(z_ - rhs.z_) > cmpLimit_) return false;
    
    return true;
}

HOST_CALLABLE void C3DVector::serialize(CSerializer& io, bool saveToStream)
{
    if(saveToStream)
    {
        io << x_;
        io << y_;
        io << z_;
        io << cmpLimit_;
    }
    else
    {
        io >> x_;
        io >> y_;
        io >> z_;
        io >> cmpLimit_;
    }
}

HOSTDEV_CALLABLE void C3DVector::copy(const C3DVector& src)
{
    x_ = src.x_;
    y_ = src.y_;
    z_ = src.z_;
    
    cmpLimit_ = src.cmpLimit_;
}

HOSTDEV_CALLABLE C3DVector C3DVector::operator=(const C3DVector& src)
{
    copy(src);
    
    return *this;
}

HOSTDEV_CALLABLE double C3DVector::operator[](int index) const
{
    if(index == 0) return x_;
    if(index == 1) return y_;
    if(index == 2) return z_;
    
    return 0.0;
}

HOSTDEV_CALLABLE bool C3DVector::isZero() const
{
    if(fabs(x_) > cmpLimit_) return false;
    if(fabs(y_) > cmpLimit_) return false;
    if(fabs(z_) > cmpLimit_) return false;
    
    return true;
}

END_CUDA_COMPATIBLE()
