#pragma once
#include "MDFFMatrices.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenCellList
{
public:
    HOSTDEV_CALLABLE CFunctorGenCellList();

public:
    void setForceFieldMatrices(CMDFFMatrices& ffMatrices);
    HOSTDEV_CALLABLE CMDFFMatrices::CCellListIndex operator()(CMDFFMatrices::CAtom& atom);
    HOST_CALLABLE void assembleCellList(mtdevice_vector<CMDFFMatrices::CCellListIndex>& devAtomCellIndices);
    HOSTDEV_CALLABLE static size_t cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t i, int maxAtomsInCell, size_t Nx, size_t Ny);
    HOSTDEV_CALLABLE static size_t cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t Nx, size_t Ny);

private:
    int* devCellList_;
    int* devCellListCount_;
    int cellCountX_;
    int cellCountY_;
    int cellCountZ_;
    float pbcWidthX_;
    float pbcWidthY_;
    float pbcWidthZ_;
    float pbcLowX_;
    float pbcLowY_;
    float pbcLowZ_;
    int maxAtomsInCell_;
};

END_CUDA_COMPATIBLE()
