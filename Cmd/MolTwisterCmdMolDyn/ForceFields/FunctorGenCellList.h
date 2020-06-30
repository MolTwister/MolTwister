#pragma once
#include "MDFFMatrices.h"
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFunctorGenCellList
{
public:
    HOSTDEV_CALLABLE CFunctorGenCellList(CMDFFMatrices::CCellList& cellList);

public:
    HOSTDEV_CALLABLE int operator()(CMDFFMatrices::CAtom& atom);
    HOSTDEV_CALLABLE static size_t cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t i, int maxAtomsInCell, size_t Nx, size_t Ny);
    HOSTDEV_CALLABLE static size_t cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t Nx, size_t Ny);

private:
    int* devCellList_;
    int* devCellListCount_;
    int cellCountX_;
    int cellCountY_;
    int cellCountZ_;
    int pbcWidthX_;
    int pbcWidthY_;
    int pbcWidthZ_;
    int pbcLowX_;
    int pbcLowY_;
    int pbcLowZ_;
    int maxAtomsInCell_;
};

END_CUDA_COMPATIBLE()
