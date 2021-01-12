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
#include "MolTwisterVariable.h"

BEGIN_CUDA_COMPATIBLE()

void CVar::serialize(CSerializer& io, bool saveToStream)
{
    int type;

    if(saveToStream)
    {
        type = (int)type_;

        io << type;
        io << name_;
    }
    else
    {
        io >> type;
        io >> name_;

        type_ = (EType)type;
    }
}

CVarAtom::CVarAtom(const CVar& src) : CVar(src)
{ 
    CVarAtom* p = (CVarAtom*)&src; 
    
    type_ = typeAtom; 
    
    if(src.getType() == typeAtom)
        atomIndex_ = p->atomIndex_; 
}

void CVarAtom::serialize(CSerializer& io, bool saveToStream)
{
    CVar::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << atomIndex_;
    }
    else
    {
        io >> atomIndex_;
    }
}

CVarBond::CVarBond(const CVar& src) : CVar(src)
{ 
    CVarBond* p = (CVarBond*)&src; 
    
    type_ = typeBond; 
    
    if(src.getType() == typeBond)
    {
        atomIndex1_ = p->atomIndex1_; 
        atomIndex2_ = p->atomIndex2_; 
    }
}

void CVarBond::serialize(CSerializer& io, bool saveToStream)
{
    CVar::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << atomIndex1_;
        io << atomIndex2_;
    }
    else
    {
        io >> atomIndex1_;
        io >> atomIndex2_;
    }
}

CVarAngle::CVarAngle(const CVar& src) : CVar(src)
{ 
    CVarAngle* p = (CVarAngle*)&src; 
    
    type_ = typeAngle; 
    
    if(src.getType() == typeAngle)
    {
        atomIndex1_ = p->atomIndex1_; 
        atomIndex2_ = p->atomIndex2_; 
        atomIndex3_ = p->atomIndex3_; 
    }
}

void CVarAngle::serialize(CSerializer& io, bool saveToStream)
{
    CVar::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << atomIndex1_;
        io << atomIndex2_;
        io << atomIndex3_;
    }
    else
    {
        io >> atomIndex1_;
        io >> atomIndex2_;
        io >> atomIndex3_;
    }
}

CVarDihedral::CVarDihedral(const CVar& src) : CVar(src)
{ 
    CVarDihedral* p = (CVarDihedral*)&src; 
    
    type_ = typeDihedral; 
    
    if(src.getType() == typeDihedral)
    {
        atomIndex1_ = p->atomIndex1_; 
        atomIndex2_ = p->atomIndex2_; 
        atomIndex3_ = p->atomIndex3_; 
        atomIndex4_ = p->atomIndex4_; 
    }
}

void CVarDihedral::serialize(CSerializer& io, bool saveToStream)
{
    CVar::serialize(io, saveToStream);

    if(saveToStream)
    {
        io << atomIndex1_;
        io << atomIndex2_;
        io << atomIndex3_;
        io << atomIndex4_;
    }
    else
    {
        io >> atomIndex1_;
        io >> atomIndex2_;
        io >> atomIndex3_;
        io >> atomIndex4_;
    }
}

END_CUDA_COMPATIBLE()
