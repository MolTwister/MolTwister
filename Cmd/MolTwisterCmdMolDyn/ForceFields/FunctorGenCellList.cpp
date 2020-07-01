#include "FunctorGenCellList.h"

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorGenCellList::CFunctorGenCellList()
{
    devCellList_ = nullptr;
    devCellListCount_ = nullptr;
}

void CFunctorGenCellList::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devCellList_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellList_[0]);
    devCellListCount_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellListCount_[0]);

    cellCountX_ = ffMatrices.cellList_. getCellCountX();
    cellCountY_ = ffMatrices.cellList_.getCellCountY();
    cellCountZ_ = ffMatrices.cellList_.getCellCountZ();

    pbcWidthX_ = ffMatrices.cellList_.getPBCWidthX();
    pbcWidthY_ = ffMatrices.cellList_.getPBCWidthY();
    pbcWidthZ_ = ffMatrices.cellList_.getPBCWidthZ();

    pbcLowX_ = ffMatrices.cellList_.getPBCLowX();
    pbcLowY_ = ffMatrices.cellList_.getPBCLowY();
    pbcLowZ_ = ffMatrices.cellList_.getPBCLowZ();

    maxAtomsInCell_ = ffMatrices.cellList_.getMaxAtomsInCell();
}

HOSTDEV_CALLABLE CMDFFMatrices::CCellListIndex CFunctorGenCellList::operator()(CMDFFMatrices::CAtom& atom)
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
    int ix = floor((r.x_ - lx) * float(Nx) / wx);
    int iy = floor((r.y_ - ly) * float(Ny) / wy);
    int iz = floor((r.z_ - lz) * float(Nz) / wz);

    ix = (ix < 0) ? 0 : ((ix >= Nx) ? Nx - 1 : ix);
    iy = (iy < 0) ? 0 : ((iy >= Ny) ? Ny - 1 : iy);
    iz = (iz < 0) ? 0 : ((iz >= Nz) ? Nz - 1 : iz);

    size_t cellIndex = cellIndexToFlatIndex((size_t)ix, (size_t)iy, (size_t)iz, Nx, Ny);
    int currentCount = devCellListCount_[cellIndex];
    if(currentCount < maxAtomsInCell_)
    {
        devCellList_[cellIndexToFlatIndex((size_t)ix, (size_t)iy, (size_t)iz, (size_t)currentCount, maxAtomsInCell_, Nx, Ny)] = (int)atom.index_;
        devCellListCount_[cellIndex]++;
    }

    return CMDFFMatrices::CCellListIndex(ix, iy, iz);
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
