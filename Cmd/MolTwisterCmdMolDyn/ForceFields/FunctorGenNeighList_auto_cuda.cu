#include "FunctorGenNeighList.h"
#include "FunctorGenCellList.h"

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorGenNeighList::CFunctorGenNeighList()
{
    devCellList_ = nullptr;
    devCellListEntryPointers_ = nullptr;
    devAtomCellIndicesRaw_ = nullptr;
    devNeighList_ = nullptr;
    devNeighListCount_ = nullptr;
    devAtomList_ = nullptr;
    numAtoms_ = 0;
    cellCountX_ = 0;
    cellCountY_ = 0;
    cellCountZ_ = 0;
    maxNeighbors_ = 0;
    rCutoff2_ = 0.0f;
}

void CFunctorGenNeighList::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devCellList_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellList_[0]);
    devCellListEntryPointers_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellListEntryPointers_[0]);
    devAtomCellIndicesRaw_ = mtraw_pointer_cast(&ffMatrices.cellList_.devAtomCellIndices_[0]);

    devNeighList_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighList_[0]);
    devNeighListCount_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighListCount_[0]);

    devAtomList_ = mtraw_pointer_cast(&ffMatrices.devAtomList_[0]);
    numAtoms_ = (int)ffMatrices.devAtomList_.size();

    cellCountX_ = ffMatrices.cellList_.getCellCountX();
    cellCountY_ = ffMatrices.cellList_.getCellCountY();
    cellCountZ_ = ffMatrices.cellList_.getCellCountZ();

    maxNeighbors_ = ffMatrices.neighList_.getMaxNeighbors();

    rCutoff2_ = ffMatrices.getRCutoff();
    rCutoff2_*= rCutoff2_;

    pbc_.rLow_.x_ = (double)ffMatrices.cellList_.getPBCLowX();
    pbc_.rLow_.y_ = (double)ffMatrices.cellList_.getPBCLowY();
    pbc_.rLow_.z_ = (double)ffMatrices.cellList_.getPBCLowZ();

    pbc_.rHigh_.x_ = double(ffMatrices.cellList_.getPBCLowX() + ffMatrices.cellList_.getPBCWidthX());
    pbc_.rHigh_.y_ = double(ffMatrices.cellList_.getPBCLowY() + ffMatrices.cellList_.getPBCWidthY());
    pbc_.rHigh_.z_ = double(ffMatrices.cellList_.getPBCLowZ() + ffMatrices.cellList_.getPBCWidthZ());
}

HOSTDEV_CALLABLE int CFunctorGenNeighList::operator()(CMDFFMatrices::CAtom& atom)
{
    if(atom.index_ >= numAtoms_) return -1;
    CMDFFMatrices::CCellListIndex cellIndex = devAtomCellIndicesRaw_[atom.index_];

    int maxCellListSize = cellCountX_ * cellCountY_ * cellCountZ_;

    int ixl = cellIndex.ix_ - 1;
    int ixh = ixl + 2;

    int iyl = cellIndex.iy_ - 1;
    int iyh = iyl + 2;

    int izl = cellIndex.iz_ - 1;
    int izh = izl + 2;

    C3DVector r0 = atom.r_;

    size_t Ix, Iy, Iz;
    for(int ix=ixl; ix<=ixh; ix++)
    {
        Ix = size_t((ix >= 0) ? ((ix < cellCountX_) ? ix : 0) : cellCountX_ - 1);
        for(int iy=iyl; iy<=iyh; iy++)
        {
            Iy = size_t((iy >= 0) ? ((iy < cellCountY_) ? iy : 0) : cellCountY_ - 1);
            for(int iz=izl; iz<=izh; iz++)
            {
                Iz = size_t((iz >= 0) ? ((iz < cellCountZ_) ? iz : 0) : cellCountZ_ - 1);
                size_t cellIndex = CFunctorGenCellList::cellIndexToFlatIndex(Ix, Iy, Iz, (size_t)cellCountX_, (size_t)cellCountY_);
                if((int)cellIndex < maxCellListSize)
                {
                    int numAtomsInCell = devCellListEntryPointers_[cellIndex].numEntries_;
                    int firstIndexInCell = devCellListEntryPointers_[cellIndex].indexFirstEntry_;
                    for(size_t i=0; i<(size_t)numAtomsInCell; i++)
                    {
                        size_t cellListEntryIndex = firstIndexInCell + i;
                        if((int)cellListEntryIndex < numAtoms_)
                        {
                            int atomIndex = devCellList_[cellListEntryIndex];
                            if(atomIndex < numAtoms_)
                            {
                                double drSqr = r0.distToAcrossPBC2(devAtomList_[atomIndex].r_, pbc_);

                                int currentNeigh = devNeighListCount_[atom.index_];
                                if(currentNeigh < maxNeighbors_)
                                {
                                    if(drSqr < (double)rCutoff2_)
                                    {
                                        if(atomIndex == atom.index_) continue;
                                        devNeighList_[neighIndexToFlatIndex(atom.index_, currentNeigh, maxNeighbors_)] = atomIndex;
                                        devNeighListCount_[atom.index_]++;
                                    }
                                }
                                else
                                {
                                    // Return the negative of the minimum number of neighbour cells required to complete the neighbour list entry
                                    return -(currentNeigh + 1);
                                }
                            }
                            else
                            {
                                return 0x00FFFFFF;
                            }
                        }
                        else
                        {
                            return 0x00FFFFFE;
                        }
                    }
                }
                else
                {
                    return 0x00FFFFFD;
                }
            }
        }
    }

    return devNeighListCount_[atom.index_];
}

HOSTDEV_CALLABLE size_t CFunctorGenNeighList::neighIndexToFlatIndex(size_t atomIndex, size_t neighIndex, int maxNeighbors)
{
    return neighIndex + maxNeighbors*atomIndex;
}

CUDA_GLOBAL void kernelCheckRequiredMaxNumNeighbours(int* devNeighListCount, int numAtoms, int* requiredMaxNeighbours)
{
    long lIdx = mtblockDim.x*mtblockIdx.x + mtthreadIdx.x;
    if(lIdx > 0) return;

    *requiredMaxNeighbours = -1;
    for(int i=0; i<numAtoms; i++)
    {
        // We found a neigh list entry that was too short
        if(devNeighListCount[i] < 0)
        {
            int minReuquiredLen = -devNeighListCount[i];
            if(minReuquiredLen > (*requiredMaxNeighbours)) *requiredMaxNeighbours = minReuquiredLen;
        }
    }
}

int CFunctorGenNeighList::checkRequiredMaxNumNeighbours(mtdevice_vector<int>& devNeighListCount)
{
    // Iterate through all lengths that was found in devNeighListCount. If negative, this means that a large enough
    // size was not found. The minimum required size is the absolute of this value. We find and return the largest
    // of these values. If we find that all neigh. list entries are large enough we return -1.
    int* devRequiredMaxNeighbours;
    mtcudaMalloc((void**)&devRequiredMaxNeighbours, sizeof(int));
    int* devNeighListCountRaw = mtraw_pointer_cast(&devNeighListCount[0]);
    kernelCheckRequiredMaxNumNeighbours CUDA_PROC(1, 1)(devNeighListCountRaw, numAtoms_, devRequiredMaxNeighbours);
    mtcudaDeviceSynchronize();

    int requiredMaxNeighbours;
    mtcudaMemcpy(&requiredMaxNeighbours, devRequiredMaxNeighbours, sizeof(int), cudaMemcpyDeviceToHost);
    mtcudaFree(devRequiredMaxNeighbours);

    return requiredMaxNeighbours;
}

END_CUDA_COMPATIBLE()
