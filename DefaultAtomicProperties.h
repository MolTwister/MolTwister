#pragma once
#include <vector>
#include <string.h>
#include <string>
#include "Utilities/3DVector.h"
#include "Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CDefaultAtomicProperties
{
public:
    class CAtomicProperty
    {
    public:
        CAtomicProperty() { ID_[0] = '\0'; }
        CAtomicProperty(const char* ID, C3DVector color, double sigma=1.0, double RCov=0.8) { strcpy(ID_, ID); color_ = color; sigma_ = sigma; RCov_ = RCov; }
        
    public:
        C3DVector color_;
        char ID_[10];
        double sigma_;
        double RCov_;
    };
    
public:
    CDefaultAtomicProperties();
    ~CDefaultAtomicProperties();
    
public:
    void getCPKColor(std::string ID, double& r, double& g, double& b) const;
    double getWDWRadius(std::string ID) const;
    double getCovalentRadius(std::string ID) const;
    int identifyAtom(std::string ID) const;
    std::string getRecognizedAtomType(std::string ID) const;

private:
    void initAtomicPropertiesArray();
    
private:
    std::vector<CAtomicProperty> atomicProperties_;
};

END_CUDA_COMPATIBLE()
