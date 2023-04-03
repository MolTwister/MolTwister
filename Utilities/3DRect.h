//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#pragma once
#include "CUDAGeneralizations.h"
#include "Serializer.h"
#include "3DVector.h"

BEGIN_CUDA_COMPATIBLE()

class C3DRect
{
public:
    HOSTDEV_CALLABLE C3DRect() {}
    HOSTDEV_CALLABLE C3DRect(C3DVector rLow, C3DVector rHigh) { rLow_ = rLow; rHigh_ = rHigh; }

public:
    HOST_CALLABLE void serialize(CSerializer& io, bool saveToStream);
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
