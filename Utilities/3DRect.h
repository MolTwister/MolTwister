#pragma once
#include "CUDAGeneralizations.h"
#include "3DVector.h"
#include <iostream>
#include <sstream>

BEGIN_CUDA_COMPATIBLE()

class C3DRect
{
public:
    HOSTDEV_CALLABLE C3DRect() {}
    HOSTDEV_CALLABLE C3DRect(C3DVector rLow, C3DVector rHigh) { rLow_ = rLow; rHigh_ = rHigh; }

public:
    HOST_CALLABLE void serialize(std::stringstream& io, bool saveToStream);
    HOSTDEV_CALLABLE double getWidthX() const { return rHigh_.x_ - rLow_.x_; }
    HOSTDEV_CALLABLE double getWidthY() const { return rHigh_.y_ - rLow_.y_; }
    HOSTDEV_CALLABLE double getWidthZ() const { return rHigh_.z_ - rLow_.z_; }
    HOSTDEV_CALLABLE double getVolume() const { return (getWidthX() * getWidthY() * getWidthZ()); }
    HOSTDEV_CALLABLE C3DVector getCenter() const;
    HOSTDEV_CALLABLE double getLargestWidth() const;
    HOSTDEV_CALLABLE void expandByFactor(double factor);
    HOSTDEV_CALLABLE void expandByLength(double length);
    HOSTDEV_CALLABLE void shrinkByLength(double length);
    HOSTDEV_CALLABLE bool isWithin(C3DVector r) const;

public:
    C3DVector rLow_;
    C3DVector rHigh_;
};

END_CUDA_COMPATIBLE()
