#include <iostream>
#include "Utilities/ASCIIUtility.h"
#include "DefaultAtomicProperties.h"

BEGIN_CUDA_COMPATIBLE()

CDefaultAtomicProperties::CDefaultAtomicProperties()
{
    atomicProperties_.reserve(200);

    initAtomicPropertiesArray();
}

CDefaultAtomicProperties::~CDefaultAtomicProperties()
{
    atomicProperties_.empty();
}

void CDefaultAtomicProperties::initAtomicPropertiesArray()
{
    // Initialize CPK colors, van der Waals radius (always place two character ID's berfore single 
    // character ID's). Most radii have been taken from [Dalton Trans., 2013, 42, 8617]. Radii that were
    // not found have been set to 1.0Å. Covalent radii are taken from [Chem. Eur. J. 2009, 15, 186].
    // Radii not found there are set to 0.5Å.
    atomicProperties_.emplace_back(CAtomicProperty("Uut", C3DVector(0.5, 0.0, 0.5), 1.00, 1.36)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Uup", C3DVector(0.5, 0.0, 0.5), 1.00, 1.62)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Uus", C3DVector(0.5, 0.0, 0.5), 1.00, 1.65)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Uuo", C3DVector(0.5, 0.0, 0.5), 1.00, 1.57)); // vdW def. 1.0
    
    atomicProperties_.emplace_back(CAtomicProperty("He", C3DVector(0.0, 0.8, 0.8), 1.43, 0.46));
    atomicProperties_.emplace_back(CAtomicProperty("Li", C3DVector(0.8, 0.0, 0.8), 2.12, 1.33));
    atomicProperties_.emplace_back(CAtomicProperty("Be", C3DVector(0.0, 0.2, 0.0), 1.98, 1.02));
    atomicProperties_.emplace_back(CAtomicProperty("Ne", C3DVector(0.0, 0.8, 0.8), 1.58, 0.67));
    atomicProperties_.emplace_back(CAtomicProperty("Na", C3DVector(0.8, 0.0, 0.8), 2.50, 1.55));
    atomicProperties_.emplace_back(CAtomicProperty("Mg", C3DVector(0.0, 0.2, 0.0), 2.51, 1.39));
    atomicProperties_.emplace_back(CAtomicProperty("Al", C3DVector(0.5, 0.0, 0.5), 2.25, 1.26));
    atomicProperties_.emplace_back(CAtomicProperty("Si", C3DVector(0.3, 0.3, 0.3), 2.19, 1.16));
    atomicProperties_.emplace_back(CAtomicProperty("Cl", C3DVector(0.0, 0.6, 0.0), 1.82, 0.99));
    atomicProperties_.emplace_back(CAtomicProperty("Ar", C3DVector(0.0, 0.8, 0.8), 1.83, 0.96));
    atomicProperties_.emplace_back(CAtomicProperty("Ca", C3DVector(0.0, 0.2, 0.0), 2.62, 1.71));
    atomicProperties_.emplace_back(CAtomicProperty("Sc", C3DVector(0.5, 0.0, 0.5), 2.58, 1.48));
    atomicProperties_.emplace_back(CAtomicProperty("Ti", C3DVector(0.3, 0.3, 0.3), 2.46, 1.36));
    atomicProperties_.emplace_back(CAtomicProperty("Cr", C3DVector(0.5, 0.0, 0.5), 2.45, 1.22));
    atomicProperties_.emplace_back(CAtomicProperty("Mn", C3DVector(0.5, 0.0, 0.5), 2.45, 1.19));
    atomicProperties_.emplace_back(CAtomicProperty("Fe", C3DVector(0.6, 0.2, 0.0), 2.44, 1.16));
    atomicProperties_.emplace_back(CAtomicProperty("Co", C3DVector(0.5, 0.0, 0.5), 2.40, 1.11));
    atomicProperties_.emplace_back(CAtomicProperty("Ni", C3DVector(0.5, 0.0, 0.5), 2.40, 1.10));
    atomicProperties_.emplace_back(CAtomicProperty("Cu", C3DVector(0.5, 0.0, 0.5), 2.38, 1.12));
    atomicProperties_.emplace_back(CAtomicProperty("Zn", C3DVector(0.5, 0.0, 0.5), 2.39, 1.18));
    atomicProperties_.emplace_back(CAtomicProperty("Ga", C3DVector(0.5, 0.0, 0.5), 2.32, 1.24));
    atomicProperties_.emplace_back(CAtomicProperty("Ge", C3DVector(0.5, 0.0, 0.5), 2.29, 1.21));
    atomicProperties_.emplace_back(CAtomicProperty("As", C3DVector(0.5, 0.0, 0.5), 1.88, 1.21));
    atomicProperties_.emplace_back(CAtomicProperty("Se", C3DVector(0.5, 0.0, 0.5), 1.82, 1.16));
    atomicProperties_.emplace_back(CAtomicProperty("Br", C3DVector(0.3, 0.0, 0.0), 1.86, 1.14));
    atomicProperties_.emplace_back(CAtomicProperty("Kr", C3DVector(0.0, 0.8, 0.8), 2.25, 1.17));
    atomicProperties_.emplace_back(CAtomicProperty("Rb", C3DVector(0.8, 0.0, 0.8), 3.21, 2.10));
    atomicProperties_.emplace_back(CAtomicProperty("Sr", C3DVector(0.0, 0.2, 0.0), 2.84, 1.85));
    atomicProperties_.emplace_back(CAtomicProperty("Zr", C3DVector(0.5, 0.0, 0.5), 2.52, 1.54));
    atomicProperties_.emplace_back(CAtomicProperty("Nb", C3DVector(0.5, 0.0, 0.5), 2.56, 1.47));
    atomicProperties_.emplace_back(CAtomicProperty("Mo", C3DVector(0.5, 0.0, 0.5), 2.45, 1.38));
    atomicProperties_.emplace_back(CAtomicProperty("Tc", C3DVector(0.5, 0.0, 0.5), 2.44, 1.28));
    atomicProperties_.emplace_back(CAtomicProperty("Ru", C3DVector(0.5, 0.0, 0.5), 2.46, 1.25));
    atomicProperties_.emplace_back(CAtomicProperty("Rh", C3DVector(0.5, 0.0, 0.5), 2.44, 1.25));
    atomicProperties_.emplace_back(CAtomicProperty("Pd", C3DVector(0.5, 0.0, 0.5), 2.15, 1.20));
    atomicProperties_.emplace_back(CAtomicProperty("Ag", C3DVector(0.5, 0.0, 0.5), 2.53, 1.28));
    atomicProperties_.emplace_back(CAtomicProperty("Cd", C3DVector(0.5, 0.0, 0.5), 2.49, 1.36));
    atomicProperties_.emplace_back(CAtomicProperty("In", C3DVector(0.5, 0.0, 0.5), 2.43, 1.42));
    atomicProperties_.emplace_back(CAtomicProperty("Sn", C3DVector(0.5, 0.0, 0.5), 2.42, 1.40));
    atomicProperties_.emplace_back(CAtomicProperty("Sb", C3DVector(0.5, 0.0, 0.5), 2.47, 1.40));
    atomicProperties_.emplace_back(CAtomicProperty("Te", C3DVector(0.5, 0.0, 0.5), 1.99, 1.36));
    atomicProperties_.emplace_back(CAtomicProperty("Xe", C3DVector(0.0, 0.8, 0.8), 2.06, 1.31));
    atomicProperties_.emplace_back(CAtomicProperty("Cs", C3DVector(0.8, 0.0, 0.8), 3.48, 2.32));
    atomicProperties_.emplace_back(CAtomicProperty("Ba", C3DVector(0.0, 0.2, 0.0), 3.03, 1.96));
    atomicProperties_.emplace_back(CAtomicProperty("Hf", C3DVector(0.5, 0.0, 0.5), 2.63, 1.52));
    atomicProperties_.emplace_back(CAtomicProperty("Ta", C3DVector(0.5, 0.0, 0.5), 2.53, 1.46));
    atomicProperties_.emplace_back(CAtomicProperty("Re", C3DVector(0.5, 0.0, 0.5), 2.49, 1.31));
    atomicProperties_.emplace_back(CAtomicProperty("Os", C3DVector(0.5, 0.0, 0.5), 2.48, 1.29));
    atomicProperties_.emplace_back(CAtomicProperty("Ir", C3DVector(0.5, 0.0, 0.5), 2.41, 1.22));
    atomicProperties_.emplace_back(CAtomicProperty("Pt", C3DVector(0.5, 0.0, 0.5), 2.29, 1.23));
    atomicProperties_.emplace_back(CAtomicProperty("Au", C3DVector(0.5, 0.0, 0.5), 2.32, 1.24));
    atomicProperties_.emplace_back(CAtomicProperty("Hg", C3DVector(0.5, 0.0, 0.5), 2.45, 1.33));
    atomicProperties_.emplace_back(CAtomicProperty("TI", C3DVector(0.5, 0.0, 0.5), 2.47, 1.44));
    atomicProperties_.emplace_back(CAtomicProperty("Pb", C3DVector(0.5, 0.0, 0.5), 2.60, 1.44));
    atomicProperties_.emplace_back(CAtomicProperty("Bi", C3DVector(0.5, 0.0, 0.5), 2.54, 1.51));
    atomicProperties_.emplace_back(CAtomicProperty("Po", C3DVector(0.5, 0.0, 0.5), 1.00, 1.45)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("At", C3DVector(0.5, 0.0, 0.5), 1.00, 1.47)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Rn", C3DVector(0.5, 0.0, 0.5), 1.00, 1.42)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Fr", C3DVector(0.8, 0.0, 0.8), 1.00, 2.23)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Ra", C3DVector(0.0, 0.2, 0.0), 1.00, 2.01)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Rf", C3DVector(0.5, 0.0, 0.5), 1.00, 1.57)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Db", C3DVector(0.5, 0.0, 0.5), 1.00, 1.49)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Sg", C3DVector(0.5, 0.0, 0.5), 1.00, 1.43)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Bh", C3DVector(0.5, 0.0, 0.5), 1.00, 1.41)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Hs", C3DVector(0.5, 0.0, 0.5), 1.00, 1.34)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Mt", C3DVector(0.5, 0.0, 0.5), 1.00, 1.29)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Ds", C3DVector(0.5, 0.0, 0.5), 1.00, 1.28)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Rg", C3DVector(0.5, 0.0, 0.5), 1.00, 1.21)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Cn", C3DVector(0.5, 0.0, 0.5), 1.00, 1.22)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("FI", C3DVector(0.5, 0.0, 0.5), 1.00, 1.43)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Lv", C3DVector(0.5, 0.0, 0.5), 1.00, 1.75)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("La", C3DVector(0.5, 0.0, 0.5), 2.98, 1.80));
    atomicProperties_.emplace_back(CAtomicProperty("Ce", C3DVector(0.5, 0.0, 0.5), 2.88, 1.63));
    atomicProperties_.emplace_back(CAtomicProperty("Pr", C3DVector(0.5, 0.0, 0.5), 2.92, 1.76));
    atomicProperties_.emplace_back(CAtomicProperty("Nd", C3DVector(0.5, 0.0, 0.5), 2.95, 1.74));
    atomicProperties_.emplace_back(CAtomicProperty("Pm", C3DVector(0.5, 0.0, 0.5), 1.00, 1.73)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Sm", C3DVector(0.5, 0.0, 0.5), 2.90, 1.72));
    atomicProperties_.emplace_back(CAtomicProperty("Eu", C3DVector(0.5, 0.0, 0.5), 2.87, 1.68));
    atomicProperties_.emplace_back(CAtomicProperty("Gd", C3DVector(0.5, 0.0, 0.5), 2.83, 1.69));
    atomicProperties_.emplace_back(CAtomicProperty("Tb", C3DVector(0.5, 0.0, 0.5), 2.79, 1.68));
    atomicProperties_.emplace_back(CAtomicProperty("Dy", C3DVector(0.5, 0.0, 0.5), 2.87, 1.67));
    atomicProperties_.emplace_back(CAtomicProperty("Ho", C3DVector(0.5, 0.0, 0.5), 2.81, 1.66));
    atomicProperties_.emplace_back(CAtomicProperty("Er", C3DVector(0.5, 0.0, 0.5), 2.83, 1.65));
    atomicProperties_.emplace_back(CAtomicProperty("Tm", C3DVector(0.5, 0.0, 0.5), 2.79, 1.64));
    atomicProperties_.emplace_back(CAtomicProperty("Yb", C3DVector(0.5, 0.0, 0.5), 2.80, 1.70));
    atomicProperties_.emplace_back(CAtomicProperty("Lu", C3DVector(0.5, 0.0, 0.5), 2.74, 1.62));
    atomicProperties_.emplace_back(CAtomicProperty("Ac", C3DVector(0.5, 0.0, 0.5), 2.80, 1.86));
    atomicProperties_.emplace_back(CAtomicProperty("Th", C3DVector(0.5, 0.0, 0.5), 2.93, 1.75));
    atomicProperties_.emplace_back(CAtomicProperty("Pa", C3DVector(0.5, 0.0, 0.5), 2.88, 1.69));
    atomicProperties_.emplace_back(CAtomicProperty("Np", C3DVector(0.5, 0.0, 0.5), 2.82, 1.71));
    atomicProperties_.emplace_back(CAtomicProperty("Pu", C3DVector(0.5, 0.0, 0.5), 2.81, 1.72));
    atomicProperties_.emplace_back(CAtomicProperty("Am", C3DVector(0.5, 0.0, 0.5), 2.83, 1.66));
    atomicProperties_.emplace_back(CAtomicProperty("Cm", C3DVector(0.5, 0.0, 0.5), 3.05, 1.66));
    atomicProperties_.emplace_back(CAtomicProperty("Bk", C3DVector(0.5, 0.0, 0.5), 3.40, 1.68));
    atomicProperties_.emplace_back(CAtomicProperty("Cf", C3DVector(0.5, 0.0, 0.5), 3.05, 1.68));
    atomicProperties_.emplace_back(CAtomicProperty("Es", C3DVector(0.5, 0.0, 0.5), 2.70, 1.65));
    atomicProperties_.emplace_back(CAtomicProperty("Fm", C3DVector(0.5, 0.0, 0.5), 1.00, 1.67)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Md", C3DVector(0.5, 0.0, 0.5), 1.00, 1.73)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("No", C3DVector(0.5, 0.0, 0.5), 1.00, 1.76)); // vdW def. 1.0
    atomicProperties_.emplace_back(CAtomicProperty("Lr", C3DVector(0.5, 0.0, 0.5), 1.00, 1.61)); // vdW def. 1.0
    
    atomicProperties_.emplace_back(CAtomicProperty("H", C3DVector(0.8, 0.8, 0.8), 1.20, 0.32));
    atomicProperties_.emplace_back(CAtomicProperty("B", C3DVector(0.9, 0.5, 0.5), 1.91, 0.85));
    atomicProperties_.emplace_back(CAtomicProperty("C", C3DVector(0.1, 0.1, 0.1), 1.77, 0.75));
    atomicProperties_.emplace_back(CAtomicProperty("N", C3DVector(0.0, 0.0, 0.8), 1.66, 0.71));
    atomicProperties_.emplace_back(CAtomicProperty("O", C3DVector(1.0, 0.0, 0.0), 1.50, 0.63));
    atomicProperties_.emplace_back(CAtomicProperty("F", C3DVector(0.0, 0.6, 0.0), 1.46, 0.64));
    atomicProperties_.emplace_back(CAtomicProperty("P", C3DVector(0.8, 0.4, 0.0), 1.90, 1.11));
    atomicProperties_.emplace_back(CAtomicProperty("S", C3DVector(0.8, 1.0, 0.0), 1.89, 1.03));
    atomicProperties_.emplace_back(CAtomicProperty("K", C3DVector(0.8, 0.0, 0.8), 2.73, 1.96));
    atomicProperties_.emplace_back(CAtomicProperty("V", C3DVector(0.5, 0.0, 0.5), 2.42, 1.34));
    atomicProperties_.emplace_back(CAtomicProperty("Y", C3DVector(0.5, 0.0, 0.5), 2.75, 1.63));
    atomicProperties_.emplace_back(CAtomicProperty("I", C3DVector(0.3, 0.0, 0.3), 2.04, 1.33));
    atomicProperties_.emplace_back(CAtomicProperty("W", C3DVector(0.5, 0.0, 0.5), 2.57, 1.37));
    atomicProperties_.emplace_back(CAtomicProperty("U", C3DVector(0.5, 0.0, 0.5), 2.71, 1.70));
}

int CDefaultAtomicProperties::identifyAtom(std::string ID) const
{
    // This function will return an ID that corresponds to an index
    // in the m_aAtomicProperties list specified in InitAtomicPropertiesArray()
    
    for(int j=0; j<(int)atomicProperties_.size(); j++)
    {
        if(CASCIIUtility::findString(atomicProperties_[j].ID_, ID) != -1)
        {
            return j;
        }
    }
    
    return -1;
}

void CDefaultAtomicProperties::getCPKColor(std::string ID, double& r, double& g, double& b) const
{
    r=0.5; g=0.0; b=0.5;

    for(int j=0; j<(int)atomicProperties_.size(); j++)
    {
        if(CASCIIUtility::findString(atomicProperties_[j].ID_, ID) != -1)
        { 
            r = atomicProperties_[j].color_.x_;
            g = atomicProperties_[j].color_.y_;
            b = atomicProperties_[j].color_.z_;
            break;
        }
    }
}

double CDefaultAtomicProperties::getWDWRadius(std::string ID) const
{
    int j = identifyAtom(ID);
    if(j != -1)
        return atomicProperties_[j].sigma_;
    
    return 1.0;
}

double CDefaultAtomicProperties::getCovalentRadius(std::string ID) const
{
    int j = identifyAtom(ID);
    if(j != -1)
        return atomicProperties_[j].RCov_;
    
    return 0.5;
}

std::string CDefaultAtomicProperties::getRecognizedAtomType(std::string ID) const
{
    int index = identifyAtom(ID);
    if((index >= 0) && (index < (int)atomicProperties_.size()))
    {
        return atomicProperties_[index].ID_;
    }
    
    return "?";
}

END_CUDA_COMPATIBLE()
