#pragma once
#include "ToolsBase.h"
#include <vector>
#include <string>
#include <map>
#include "../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CMolecularSystemTools : public CToolsBase
{
public:
    CMolecularSystemTools() = delete;
    CMolecularSystemTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    void removeDuplicateBondIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>* auxList=nullptr) const;
    void removeDuplicateAngleIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>* auxList=nullptr) const;
    void removeDuplicateDihIndices(std::vector<int>& atoms1, std::vector<int>& atoms2, std::vector<int>& atoms3, std::vector<int>& atoms4, std::vector<int>* auxList=nullptr) const;
    void genMolIndices(const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molIndices) const;
    void getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule) const;

private:
    void getMoleculeConnectedToIndex(int index, const std::vector<std::vector<int>>& bondDestIndices, std::vector<int>& molecule, std::map<int, int>& visitationTracker) const;
};

END_CUDA_COMPATIBLE()
