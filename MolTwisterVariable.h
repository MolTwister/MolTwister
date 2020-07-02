#pragma once
#include <iostream>
#include <sstream>
#include "MolTwisterAtom.h"

class CVar
{
public:
    enum EType { typeEmpty=0, typeAtom=1, typeBond=2, typeAngle=3, typeDihedral=4 };
    
    
public:
    CVar() { type_ = typeEmpty; name_[0] = '\0'; }
    CVar(const CVar& src) { type_ = src.type_; name_ = src.name_; }
    
public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::shared_ptr<CVar> createCopy() const = 0;
    virtual std::shared_ptr<CVar> createCopy(const CVar& src) const = 0;
    EType getType() const { return type_; }
    void setName(const char* name) { name_ = name; }
    void getName(std::string& name) const { name = name_; }
    
protected:
    EType type_;
    std::string name_;
};

class CVarAtom : public CVar
{
public:
    CVarAtom() : CVar() { type_ = typeAtom; atomIndex_ = -1; }
    CVarAtom(const CVar& src);

public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::shared_ptr<CVar> createCopy() const { auto ret = std::shared_ptr<CVar>(new CVarAtom); *(CVarAtom*)(ret.get()) = *this; return ret; }
    virtual std::shared_ptr<CVar> createCopy(const CVar& src) const { auto ret = std::shared_ptr<CVar>(new CVarAtom(src)); *(CVarAtom*)(ret.get()) = *this; return ret; }

public:
    int atomIndex_;
};

class CVarBond : public CVar
{
public:
    CVarBond() : CVar() { type_ = typeBond; atomIndex1_ = -1; atomIndex2_ = -1; }
    CVarBond(const CVar& src);

public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::shared_ptr<CVar> createCopy() const { auto ret = std::shared_ptr<CVar>(new CVarBond); *(CVarBond*)(ret.get()) = *this; return ret; }
    virtual std::shared_ptr<CVar> createCopy(const CVar& src) const { auto ret = std::shared_ptr<CVar>(new CVarBond(src)); *(CVarBond*)(ret.get()) = *this; return ret; }

public:
    int atomIndex1_;
    int atomIndex2_;
};

class CVarAngle : public CVar
{
public:
    CVarAngle() : CVar() { type_ = typeAngle; atomIndex1_ = -1; atomIndex2_ = -1; atomIndex3_ = -1; }
    CVarAngle(const CVar& src);

public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::shared_ptr<CVar> createCopy() const { auto ret = std::shared_ptr<CVar>(new CVarAngle); *(CVarAngle*)(ret.get()) = *this; return ret; }
    virtual std::shared_ptr<CVar> createCopy(const CVar& src) const { auto ret = std::shared_ptr<CVar>(new CVarAngle(src)); *(CVarAngle*)(ret.get()) = *this; return ret; }

public:
    int atomIndex1_;
    int atomIndex2_;
    int atomIndex3_;
};

class CVarDihedral : public CVar
{
public:
    CVarDihedral() : CVar() { type_ = typeDihedral; atomIndex1_ = -1; atomIndex2_ = -1; atomIndex3_ = -1; atomIndex4_ = -1; }
    CVarDihedral(const CVar& src);

public:
    virtual void serialize(std::stringstream& io, bool saveToStream);
    virtual std::shared_ptr<CVar> createCopy() const { auto ret = std::shared_ptr<CVar>(new CVarDihedral); *(CVarDihedral*)(ret.get()) = *this; return ret; }
    virtual std::shared_ptr<CVar> createCopy(const CVar& src) const { auto ret = std::shared_ptr<CVar>(new CVarDihedral(src)); *(CVarDihedral*)(ret.get()) = *this; return ret; }

public:
    int atomIndex1_;
    int atomIndex2_;
    int atomIndex3_;
    int atomIndex4_;
};
