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

#include "3DRect.h"

BEGIN_CUDA_COMPATIBLE()

HOST_CALLABLE void C3DRect::serialize(CSerializer& io, bool saveToStream)
{
    if(saveToStream)
    {
        rLow_.serialize(io, saveToStream);
        rHigh_.serialize(io, saveToStream);
    }
    else
    {
        rLow_.serialize(io, saveToStream);
        rHigh_.serialize(io, saveToStream);
    }
}

HOSTDEV_CALLABLE C3DVector C3DRect::getCenter() const
{
    C3DVector center((rLow_.x_ + rHigh_.x_) / 2.0,
                     (rLow_.y_ + rHigh_.y_) / 2.0,
                     (rLow_.z_ + rHigh_.z_) / 2.0);

    return center;
}

HOSTDEV_CALLABLE double C3DRect::getLargestWidth() const
{
    double largest = getWidthX();

    if(getWidthY() > largest) largest = getWidthY();
    if(getWidthZ() > largest) largest = getWidthZ();

    return largest;
}

HOSTDEV_CALLABLE void C3DRect::expandByFactor(double factor)
{
    double expX = getWidthX() * factor;
    double expY = getWidthY() * factor;
    double expZ = getWidthZ() * factor;
    double expXHalf = expX / 2.0;
    double expYHalf = expY / 2.0;
    double expZHalf = expZ / 2.0;

    rLow_.x_-= expXHalf;
    rLow_.y_-= expYHalf;
    rLow_.z_-= expZHalf;

    rHigh_.x_+= expXHalf;
    rHigh_.y_+= expYHalf;
    rHigh_.z_+= expZHalf;
}

HOSTDEV_CALLABLE void C3DRect::expandByLength(double length)
{
    rLow_.x_-= length;
    rLow_.y_-= length;
    rLow_.z_-= length;

    rHigh_.x_+= length;
    rHigh_.y_+= length;
    rHigh_.z_+= length;
}

HOSTDEV_CALLABLE void C3DRect::shrinkByLength(double length)
{
    rLow_.x_+= length;
    rLow_.y_+= length;
    rLow_.z_+= length;

    rHigh_.x_-= length;
    rHigh_.y_-= length;
    rHigh_.z_-= length;
}

HOSTDEV_CALLABLE bool C3DRect::isWithin(C3DVector r) const
{
    if((r.x_ >= rLow_.x_) && ((r.x_ <= rHigh_.x_)))
    {
        if((r.y_ >= rLow_.y_) && ((r.y_ <= rHigh_.y_)))
        {
            if((r.z_ >= rLow_.z_) && ((r.z_ <= rHigh_.z_)))
            {
                return true;
            }
        }
    }

    return false;
}

END_CUDA_COMPATIBLE()
