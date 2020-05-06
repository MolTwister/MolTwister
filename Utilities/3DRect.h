#pragma once
#include "3DVector.h"

class C3DRect
{
public:
    C3DRect() {}
    C3DRect(C3DVector rLow, C3DVector rHigh) { rLow_ = rLow; rHigh_ = rHigh; }

public:
    double getWidthX() const { return rHigh_.x_ - rLow_.x_; }
    double getWidthY() const { return rHigh_.y_ - rLow_.y_; }
    double getWidthZ() const { return rHigh_.z_ - rLow_.z_; }
    double getVolume() const { return (getWidthX() * getWidthY() * getWidthZ()); }
    C3DVector getCenter() const;
    double getLargestWidth() const;
    void expandByFactor(double factor);
    void expandByLength(double length);
    void shrinkByLength(double length);
    bool isWithin(C3DVector r) const;

public:
    C3DVector rLow_;
    C3DVector rHigh_;
};
