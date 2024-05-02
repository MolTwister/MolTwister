//
// Copyright (C) 2021 Richard Olsen.
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

#include <iostream>
#include "MolTwisterGLObject.h"

BEGIN_CUDA_COMPATIBLE()

void CGLObject::serialize(CSerializer& io, bool saveToStream)
{
    int type;

    if(saveToStream)
    {
        type = (int)type_;
        io << type;
    }
    else
    {
        io >> type;
        type_ = (EType)type;
    }
}

CGLObjectLine::CGLObjectLine(C3DVector p1, C3DVector p2, float r, float g, float b) : CGLObject()
{ 
    type_ = objLine;
    p1_ = p1; 
    p2_ = p2;
    
    color_[0] = r;
    color_[1] = g;
    color_[2] = b;
}

CGLObjectLine::CGLObjectLine(const CGLObject& src) : CGLObject()
{
    CGLObjectLine*  p;
    
    if(src.getType() == objLine)
    {
        p = (CGLObjectLine*)&src;
        
        type_ = objLine;
        *this = *p;
    }
}

void CGLObjectLine::serialize(CSerializer& io, bool saveToStream)
{
    CGLObject::serialize(io, saveToStream);

    if(saveToStream)
    {
        p1_.serialize(io, saveToStream);
        p2_.serialize(io, saveToStream);
        io << color_[0];
        io << color_[1];
        io << color_[2];
    }
    else
    {
        p1_.serialize(io, saveToStream);
        p2_.serialize(io, saveToStream);
        io >> color_[0];
        io >> color_[1];
        io >> color_[2];
    }
}

END_CUDA_COMPATIBLE()
