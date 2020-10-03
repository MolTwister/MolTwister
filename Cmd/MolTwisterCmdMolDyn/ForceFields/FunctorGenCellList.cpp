#include "FunctorGenCellList.h"

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorGenCellList::CFunctorGenCellList()
{
    devCellList_ = nullptr;
    devCellListEntryPointers_ = nullptr;
}

void CFunctorGenCellList::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devCellList_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellList_[0]);
    devCellListEntryPointers_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellListEntryPointers_[0]);

    cellCountX_ = ffMatrices.cellList_. getCellCountX();
    cellCountY_ = ffMatrices.cellList_.getCellCountY();
    cellCountZ_ = ffMatrices.cellList_.getCellCountZ();

    pbcWidthX_ = ffMatrices.cellList_.getPBCWidthX();
    pbcWidthY_ = ffMatrices.cellList_.getPBCWidthY();
    pbcWidthZ_ = ffMatrices.cellList_.getPBCWidthZ();

    pbcLowX_ = ffMatrices.cellList_.getPBCLowX();
    pbcLowY_ = ffMatrices.cellList_.getPBCLowY();
    pbcLowZ_ = ffMatrices.cellList_.getPBCLowZ();
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
    int ix = (int)floor((r.x_ - double(lx)) * double(Nx) / double(wx));
    int iy = (int)floor((r.y_ - double(ly)) * double(Ny) / double(wy));
    int iz = (int)floor((r.z_ - double(lz)) * double(Nz) / double(wz));

    ix = (ix < 0) ? 0 : ((ix >= Nx) ? Nx - 1 : ix);
    iy = (iy < 0) ? 0 : ((iy >= Ny) ? Ny - 1 : iy);
    iz = (iz < 0) ? 0 : ((iz >= Nz) ? Nz - 1 : iz);

    return CMDFFMatrices::CCellListIndex(ix, iy, iz);
}

CUDA_GLOBAL void kernelAssembleList(CMDFFMatrices::CCellListIndex* devAtomCellIndices, int numAtomCellIndices, int* devCellList, CMDFFMatrices::CListPointer* devCellListEntryPointers, int Nx, int Ny, int Nz, int* err)
{
    long lIdx = mtblockDim.x*mtblockIdx.x + mtthreadIdx.x;
    if(lIdx > 0) return;

    // Iterate through all atoms once to generate the list of entry pointers
    for(int i=0; i<numAtomCellIndices; i++)
    {
        CMDFFMatrices::CCellListIndex index3D = devAtomCellIndices[i];
        size_t ix = (size_t)index3D.ix_;
        size_t iy = (size_t)index3D.iy_;
        size_t iz = (size_t)index3D.iz_;

        size_t cellIndex = CFunctorGenCellList::cellIndexToFlatIndex(ix, iy, iz, Nx, Ny);
        devCellListEntryPointers[cellIndex].numEntries_++;
    }

    // Now place the starting indices, based on the counts that were found and reset the counters
    int totalCellCount = Nx * Ny * Nz;
    int nextStartingIndex = 0;
    for(int i=0; i<totalCellCount; i++)
    {
        if(devCellListEntryPointers[i].numEntries_ != 0)
        {
            devCellListEntryPointers[i].indexFirstEntry_ = nextStartingIndex;
            nextStartingIndex+= devCellListEntryPointers[i].numEntries_;
            devCellListEntryPointers[i].numEntries_ = 0;
        }
    }

    // Use the pre-filled list of entry pointers and place all the cell list entries
    bool cellListEntryOverflow = false;
    bool cellListCountOverflow = false;
    for(int i=0; i<numAtomCellIndices; i++)
    {
        CMDFFMatrices::CCellListIndex index3D = devAtomCellIndices[i];
        size_t ix = (size_t)index3D.ix_;
        size_t iy = (size_t)index3D.iy_;
        size_t iz = (size_t)index3D.iz_;

        size_t cellIndex = CFunctorGenCellList::cellIndexToFlatIndex(ix, iy, iz, Nx, Ny);
        if((int)cellIndex < totalCellCount)
        {
            CMDFFMatrices::CListPointer listPointer = devCellListEntryPointers[cellIndex];
            if(listPointer.indexFirstEntry_ != -1)
            {
                int cellListEntryIndex = listPointer.indexFirstEntry_ + listPointer.numEntries_;
                if(cellListEntryIndex < numAtomCellIndices)
                {
                    devCellList[cellListEntryIndex] = i;
                    devCellListEntryPointers[cellIndex].numEntries_++;
                }
                else
                {
                    cellListEntryOverflow = true;
                }
            }
        }
        else
        {
            cellListCountOverflow = true;
        }
    }

    if(cellListCountOverflow) *err = 1;
    else if(cellListEntryOverflow) *err = 2;
    else *err = 0;
}

HOST_CALLABLE void CFunctorGenCellList::assembleCellList(mtdevice_vector<CMDFFMatrices::CCellListIndex>& devAtomCellIndices)
{
    int* devErr;
    mtcudaMalloc((void**)&devErr, sizeof(int));
    CMDFFMatrices::CCellListIndex* devAtomCellIndicesRaw = mtraw_pointer_cast(&devAtomCellIndices[0]);
    kernelAssembleList CUDA_PROC(1, 1)(devAtomCellIndicesRaw, (int)devAtomCellIndices.size(), devCellList_, devCellListEntryPointers_, cellCountX_, cellCountY_, cellCountZ_, devErr);
    mtcudaDeviceSynchronize();

    int err;
    mtcudaMemcpy(&err, devErr, sizeof(int), cudaMemcpyDeviceToHost);
    mtcudaFree(devErr);

    if(err == 1)
    {
        printf("Warning: Cell list count overflow occurred!\r\n");
    }
    if(err == 2)
    {
        printf("Warning: Cell list entry overflow occurred!\r\n");
    }
}

HOSTDEV_CALLABLE size_t CFunctorGenCellList::cellIndexToFlatIndex(size_t ix, size_t iy, size_t iz, size_t Nx, size_t Ny)
{
    return ix + Nx*(iy + iz*Ny);
}

END_CUDA_COMPATIBLE()
