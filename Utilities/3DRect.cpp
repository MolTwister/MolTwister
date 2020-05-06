#include "3DRect.h"

C3DVector C3DRect::getCenter() const
{
    C3DVector   center((rLow_.x_ + rHigh_.x_) / 2.0,
                       (rLow_.y_ + rHigh_.y_) / 2.0,
                       (rLow_.z_ + rHigh_.z_) / 2.0);

    return center;
}

double C3DRect::getLargestWidth() const
{
    double  largest = getWidthX();

    if(getWidthY() > largest) largest = getWidthY();
    if(getWidthZ() > largest) largest = getWidthZ();

    return largest;
}

void C3DRect::expandByFactor(double factor)
{
    double  expX = getWidthX() * factor;
    double  expY = getWidthY() * factor;
    double  expZ = getWidthZ() * factor;
    double  expXHalf = expX / 2.0;
    double  expYHalf = expY / 2.0;
    double  expZHalf = expZ / 2.0;

    rLow_.x_-= expXHalf;
    rLow_.y_-= expYHalf;
    rLow_.z_-= expZHalf;

    rHigh_.x_+= expXHalf;
    rHigh_.y_+= expYHalf;
    rHigh_.z_+= expZHalf;
}

void C3DRect::expandByLength(double length)
{
    rLow_.x_-= length;
    rLow_.y_-= length;
    rLow_.z_-= length;

    rHigh_.x_+= length;
    rHigh_.y_+= length;
    rHigh_.z_+= length;
}

void C3DRect::shrinkByLength(double length)
{
    rLow_.x_+= length;
    rLow_.y_+= length;
    rLow_.z_+= length;

    rHigh_.x_-= length;
    rHigh_.y_-= length;
    rHigh_.z_-= length;
}

bool C3DRect::isWithin(C3DVector r) const
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
