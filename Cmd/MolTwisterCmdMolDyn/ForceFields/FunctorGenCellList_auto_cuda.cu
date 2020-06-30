#include "FunctorGenCellList.h"

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorGenCellList::CFunctorGenCellList(CMDFFMatrices::CCellList& cellList)
{
    devCellList_ = cellList.getCellListRaw();
    devCellListCount_ = cellList.getCellListCountRaw();

    cellCountX_ = cellList.getCellCountX();
    cellCountY_ = cellList.getCellCountY();
    cellCountZ_ = cellList.getCellCountZ();

    pbcWidthX_ = cellList.getPBCWidthX();
    pbcWidthY_ = cellList.getPBCWidthY();
    pbcWidthZ_ = cellList.getPBCWidthZ();

    pbcLowX_ = cellList.getPBCLowX();
    pbcLowY_ = cellList.getPBCLowY();
    pbcLowZ_ = cellList.getPBCLowZ();

    maxAtomsInCell_ = cellList.getMaxAtomsInCell();
}

HOSTDEV_CALLABLE int CFunctorGenCellList::operator()(CMDFFMatrices::CAtom& atom)
{
    float wx = pbcWidthX_;
    float wy = pbcWidthY_;
    float wz = pbcWidthZ_;

    float lx = pbcLowX_;
    float ly = pbcLowY_;
    float lz = pbcLowZ_;

    int Nx = cellCountX_;
    int Ny = cellCountY_;
    int Nz = cellCountZ_;

    C3DVector r = atom.r_;
    size_t ix = (size_t)floor((r.x_ - lx) * float(Nx) / wx);
    size_t iy = (size_t)floor((r.y_ - ly) * float(Ny) / wy);
    size_t iz = (size_t)floor((r.z_ - lz) * float(Nz) / wz);

    size_t cellIndex = cellIndexToFlatIndex(ix, iy, iz, Nx, Ny);
    size_t currentCount = devCellListCount_[cellIndex];
    if(int(currentCount) < maxAtomsInCell_)
    {
        devCellList_[cellIndexToFlatIndex(ix, iy, iz, currentCount, maxAtomsInCell_, Nx, Ny)] = (int)atom.index_;
        devCellListCount_[cellIndex]++;
    }

    return cellIndex;
}


HOSTDEV_CALLABLE size_t CFunctorGenCellList::cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t i, int maxAtomsInCell, size_t Nx, size_t Ny)
{
    return size_t(i + maxAtomsInCell*(ix + Nx*(iy + iz*Ny)));
}

HOSTDEV_CALLABLE size_t CFunctorGenCellList::cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t Nx, size_t Ny)
{
    return ix + Nx*(iy + iz*Ny);
}

END_CUDA_COMPATIBLE()
