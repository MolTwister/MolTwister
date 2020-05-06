#pragma once

class C3DRect;

class C3DVector
{
public:
    C3DVector() { x_=0.0; y_=0.0; z_=0.0; cmpLimit_=1E-16; }
    C3DVector(double x, double y, double z) { x_=x; y_=y; z_=z; cmpLimit_=1E-16; }
    C3DVector(const C3DVector& src) { copy(src); }
    
public:
    void copy(const C3DVector& src);
    void set(double x, double y, double z) { x_=x; y_=y; z_=z; }
    void set(const double* r) { x_=r[0]; y_=r[1]; z_=r[2]; }
    void get(double& x, double& y, double& z) const { x=x_; y=y_; z=z_; }
    void get(double* r) const { r[0]=x_; r[1]=y_; r[2]=z_; }
    void print() const;
    double norm() const;
    double norm2() const;
    double distTo(const C3DVector& v) const ;
    double distToAcrossPBC(const C3DVector& v, const C3DRect& pbc) const;
    double distToAcrossPBC2(const C3DVector& v, const C3DRect& pbc) const;
    double cosAngle() const;
    double cosAngle(const C3DVector& v2) const;
    double angle() const;
    double angle(const C3DVector& v2) const;
    void moveToSameSideOfPBCAsThis(C3DVector& v, const C3DRect& pbc) const;
    double cosAngleAcrossPBC(C3DVector v1, C3DVector v3, const C3DRect& pbc) const;
    double angleAcrossPBC(const C3DVector& v1, const C3DVector& v3, const C3DRect& pbc) const;
    C3DVector unit() const;
    void normalize();
    C3DVector cross(const C3DVector& rhs) const;
    C3DVector operator+(const C3DVector& rhs) const;
    C3DVector operator+=(const C3DVector& rhs);
    C3DVector operator-(const C3DVector& rhs) const;
    C3DVector operator-=(const C3DVector& rhs);
    double operator*(const C3DVector& rhs) const;
    C3DVector operator*(double a) const;
    C3DVector operator*=(double a);
    bool operator==(const C3DVector& rhs) const;
    C3DVector operator=(const C3DVector& src);
    double operator[](int index) const;
    bool isZero() const;
    
public:
    double x_;
    double y_;
    double z_;
    
private:
    double cmpLimit_;
};
