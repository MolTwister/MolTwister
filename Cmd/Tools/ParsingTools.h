#pragma once
#include "ToolsBase.h"
#include <string>

class CParsingTools : public CToolsBase
{
public:
    CParsingTools() = delete;
    CParsingTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    bool retrieve1to4BondSepCoeffs(std::vector<std::string> arguments, size_t& arg, double* a1to4BondSepCoeffs) const;
    std::pair<bool, std::string> retrieveDihedralAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, CAtom** atom4Ptr, int* atomIndices) const;
    std::pair<bool, std::string> retrieveAngleAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, int* atomIndices) const;
    std::pair<bool, std::string> retrieveBondAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, int* atomIndices) const;
};
