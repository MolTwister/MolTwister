#pragma once
#include "../../../Utilities/CUDAGeneralizations.h"
#include "MDFFMatrices.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenNeighList
{
public:
    HOSTDEV_CALLABLE CFunctorGenNeighList();

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE int operator()(CMDFFMatrices::CAtom& atom);
    HOSTDEV_CALLABLE static size_t neighIndexToFlatIndex(size_t atomIndex, size_t neighIndex, int maxNeighbors);

private:
    int* devCellList_;
    int* devCellListCount_;
    CMDFFMatrices::CCellListIndex* devAtomCellIndicesRaw_;
    int* devNeighList_;
    int* devNeighListCount_;
    int cellCountX_;
    int cellCountY_;
    int cellCountZ_;
    int maxNeighbors_;
    int maxAtomsInCell_;
};

END_CUDA_COMPATIBLE()
