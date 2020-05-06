#pragma once
#include "MolTwisterCmd.h"

class CCmdGauss9 : public CCmd
{
public:
    CCmdGauss9() = delete;
    CCmdGauss9(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "gauss9"; }
    std::string getTopLevHelpString() const { return std::string("Operations related to the Gaussian9 SW package"); }
    std::string getHelpString() const;

protected:
    virtual void onAddKeywords();
    
private:
    void parseDihedralrotCommand(std::string commandLine, int& arg);
    void parseAnglerotCommand(std::string commandLine, int& arg);
    void parseBondstretchCommand(std::string commandLine, int& arg);
    void printGaussianSPScript(std::string baseSetSpec, double charge, int spinMul, std::string heading) const;
    void printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                   int atInd1, int atInd2, int atInd3, int atInd4, double startAngle, double endAngle, double step) const;
    void printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                   int atInd1, int atInd2, int atInd3, double startAngle, double endAngle, double step) const;
    void printGaussianModRedScript(std::string baseSetSpec, double charge, int spinMul, std::string heading,
                                   int atInd1, int atInd2, double startDist, double endDist, double step) const;
    void parseGenxyzfrominputCommand(std::string commandLine, int& arg);
    void parseGenoptefromoutputCommand(std::string commandLine, int& arg);
    void parseGenoptxyzfromoutputCommand(std::string commandLine, int& arg);
    void parseGenxyzfromoutputCommand(std::string commandLine, int& arg);
    void printXYZ(const std::vector<std::string>& IDs, const std::vector<C3DVector>& atoms) const;
};
