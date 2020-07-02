#pragma once
#include <iostream>
#include <sstream>
#include "CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class C3DRect;

class C3DVector
{
public:
    HOSTDEV_CALLABLE C3DVector() { x_=0.0; y_=0.0; z_=0.0; cmpLimit_=1E-16; }
    HOSTDEV_CALLABLE C3DVector(double x, double y, double z) { x_=x; y_=y; z_=z; cmpLimit_=1E-16; }
    HOSTDEV_CALLABLE C3DVector(const C3DVector& src) { copy(src); }
    
public:
    HOST_CALLABLE void serialize(std::stringstream& io, bool saveToStream);
    HOSTDEV_CALLABLE void copy(const C3DVector& src);
    HOSTDEV_CALLABLE void set(double x, double y, double z) { x_=x; y_=y; z_=z; }
    HOSTDEV_CALLABLE void set(const double* r) { x_=r[0]; y_=r[1]; z_=r[2]; }
    HOSTDEV_CALLABLE void get(double& x, double& y, double& z) const { x=x_; y=y_; z=z_; }
    HOSTDEV_CALLABLE void get(double* r) const { r[0]=x_; r[1]=y_; r[2]=z_; }
    void print() const;
    HOSTDEV_CALLABLE double norm() const;
    HOSTDEV_CALLABLE double norm2() const;
    HOSTDEV_CALLABLE double distTo(const C3DVector& v) const ;
    HOSTDEV_CALLABLE double distToAcrossPBC(const C3DVector& v, const C3DRect& pbc) const;
    HOSTDEV_CALLABLE double distToAcrossPBC2(const C3DVector& v, const C3DRect& pbc) const;
    HOSTDEV_CALLABLE double cosAngle() const;
    HOSTDEV_CALLABLE double cosAngle(const C3DVector& v2) const;
    HOSTDEV_CALLABLE double angle() const;
    HOSTDEV_CALLABLE double angle(const C3DVector& v2) const;
    HOSTDEV_CALLABLE void moveToSameSideOfPBCAsThis(C3DVector& v, const C3DRect& pbc) const;
    HOSTDEV_CALLABLE double cosAngleAcrossPBC(C3DVector v1, C3DVector v3, const C3DRect& pbc) const;
    HOSTDEV_CALLABLE double angleAcrossPBC(const C3DVector& v1, const C3DVector& v3, const C3DRect& pbc) const;
    HOSTDEV_CALLABLE C3DVector unit() const;
    HOSTDEV_CALLABLE void normalize();
    HOSTDEV_CALLABLE C3DVector cross(const C3DVector& rhs) const;
    HOSTDEV_CALLABLE C3DVector operator+(const C3DVector& rhs) const;
    HOSTDEV_CALLABLE C3DVector operator+=(const C3DVector& rhs);
    HOSTDEV_CALLABLE C3DVector operator-(const C3DVector& rhs) const;
    HOSTDEV_CALLABLE C3DVector operator-=(const C3DVector& rhs);
    HOSTDEV_CALLABLE double operator*(const C3DVector& rhs) const;
    HOSTDEV_CALLABLE C3DVector operator*(double a) const;
    HOSTDEV_CALLABLE C3DVector operator*=(double a);
    HOSTDEV_CALLABLE bool operator==(const C3DVector& rhs) const;
    HOSTDEV_CALLABLE C3DVector operator=(const C3DVector& src);
    HOSTDEV_CALLABLE double operator[](int index) const;
    HOSTDEV_CALLABLE bool isZero() const;
    
public:
    double x_;
    double y_;
    double z_;
    
private:
    double cmpLimit_;
};

END_CUDA_COMPATIBLE()
