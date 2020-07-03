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

    return CMDFFMatrices::CCellListIndex(ix, iy, iz);
}

CUDA_GLOBAL void kernelAssembleList(CMDFFMatrices::CCellListIndex* devAtomCellIndices, int numAtomCellIndices, int* devCellList, int* devCellListCount, int Nx, int Ny, int Nz, int maxAtomsInCell, int* err)
{
    long lIdx = mtblockDim.x*mtblockIdx.x + mtthreadIdx.x;
    if(lIdx > 0) return;

    int maxCellListCountSize = Nx * Ny * Nz;
    int maxCellListSize = maxCellListCountSize * maxAtomsInCell;

    bool cellListOverflow = false;
    bool cellListEntryOverflow = false;
    bool cellListCountOverflow = false;
    for(int i=0; i<numAtomCellIndices; i++)
    {
        CMDFFMatrices::CCellListIndex index3D = devAtomCellIndices[i];
        size_t ix = (size_t)index3D.ix_;
        size_t iy = (size_t)index3D.iy_;
        size_t iz = (size_t)index3D.iz_;

        size_t cellIndex = CFunctorGenCellList::cellIndexToFlatIndex(ix, iy, iz, Nx, Ny);
        if((int)cellIndex < maxCellListCountSize)
        {
            int currentCount = devCellListCount[cellIndex];
            if(currentCount < maxAtomsInCell)
            {
                int cellListEntryIndex = CFunctorGenCellList::cellIndexToFlatIndex(ix, iy, iz, (size_t)currentCount, maxAtomsInCell, Nx, Ny);
                if(cellListEntryIndex < maxCellListSize)
                {
                    devCellList[cellListEntryIndex] = i;
                    devCellListCount[cellIndex]++;
                }
                else
                {
                    cellListEntryOverflow = true;
                }
            }
            else
            {
                cellListOverflow = true;
            }
        }
        else
        {
            cellListCountOverflow = true;
        }
    }

    if(cellListOverflow) *err = 1;
    else if(cellListCountOverflow) *err = 2;
    else if(cellListEntryOverflow) *err = 3;
    else *err = 0;
}

HOST_CALLABLE void CFunctorGenCellList::assembleCellList(mtdevice_vector<CMDFFMatrices::CCellListIndex>& devAtomCellIndices)
{
    int* devErr;
    mtcudaMalloc((void**)&devErr, sizeof(int));
    CMDFFMatrices::CCellListIndex* devAtomCellIndicesRaw = mtraw_pointer_cast(&devAtomCellIndices[0]);
    kernelAssembleList CUDA_PROC(1, 1)(devAtomCellIndicesRaw, (int)devAtomCellIndices.size(), devCellList_, devCellListCount_, cellCountX_, cellCountY_, cellCountZ_, maxAtomsInCell_, devErr);
    mtcudaDeviceSynchronize();

    int err;
    mtcudaMemcpy(&err, devErr, sizeof(int), cudaMemcpyDeviceToHost);
    mtcudaFree(devErr);

    if(err == 1)
    {
        printf("Warning: Cell list overflow occurred!\r\n");
    }
    if(err == 2)
    {
        printf("Warning: Cell list count overflow occurred!\r\n");
    }
    if(err == 3)
    {
        printf("Warning: Cell list entry overflow occurred!\r\n");
    }
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
