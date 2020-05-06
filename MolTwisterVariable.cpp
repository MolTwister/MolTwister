#include <iostream>
#include "MolTwisterVariable.h"

CVarAtom::CVarAtom(const CVar& src) : CVar(src)
{ 
    CVarAtom* p = (CVarAtom*)&src; 
    
    type_ = typeAtom; 
    
    if(src.getType() == typeAtom)
        atomIndex_ = p->atomIndex_; 
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
