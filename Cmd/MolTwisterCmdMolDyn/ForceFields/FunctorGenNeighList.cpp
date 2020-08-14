#include "FunctorGenNeighList.h"
#include "FunctorGenCellList.h"

BEGIN_CUDA_COMPATIBLE()

HOSTDEV_CALLABLE CFunctorGenNeighList::CFunctorGenNeighList()
{
    devCellList_ = nullptr;
    devCellListCount_ = nullptr;
    devAtomCellIndicesRaw_ = nullptr;
    devNeighList_ = nullptr;
    devNeighListCount_ = nullptr;
    devAtomList_ = nullptr;
    numAtoms_ = 0;
    cellCountX_ = 0;
    cellCountY_ = 0;
    cellCountZ_ = 0;
    maxNeighbors_ = 0;
    maxAtomsInCell_ = 0;
    rCutoff2_ = 0.0f;
}

void CFunctorGenNeighList::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devCellList_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellList_[0]);
    devCellListCount_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellListCount_[0]);
    devAtomCellIndicesRaw_ = mtraw_pointer_cast(&ffMatrices.cellList_.devAtomCellIndices_[0]);

    devNeighList_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighList_[0]);
    devNeighListCount_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighListCount_[0]);

    devAtomList_ = mtraw_pointer_cast(&ffMatrices.devAtomList_[0]);
    numAtoms_ = (int)ffMatrices.devAtomList_.size();

    cellCountX_ = ffMatrices.cellList_.getCellCountX();
    cellCountY_ = ffMatrices.cellList_.getCellCountY();
    cellCountZ_ = ffMatrices.cellList_.getCellCountZ();

    maxNeighbors_ = ffMatrices.neighList_.getMaxNeighbors();
    maxAtomsInCell_ = ffMatrices.cellList_.getMaxAtomsInCell();

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

    int maxCellListCountSize = cellCountX_ * cellCountY_ * cellCountZ_;
    int maxCellListSize = maxCellListCountSize * maxAtomsInCell_;

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
        Ix = size_t((ix >= 0) ? ((ix < ixh) ? ix : 0) : cellCountX_ - 1);
        for(int iy=iyl; iy<=iyh; iy++)
        {
            Iy = size_t((iy >= 0) ? ((iy < iyh) ? iy : 0) : cellCountY_ - 1);
            for(int iz=izl; iz<=izh; iz++)
            {
                Iz = size_t((iz >= 0) ? ((iz < izh) ? iz : 0) : cellCountZ_ - 1);
                size_t cellIndex = CFunctorGenCellList::cellIndexToFlatIndex(Ix, Iy, Iz, (size_t)cellCountX_, (size_t)cellCountY_);
                if((int)cellIndex < maxCellListCountSize)
                {
                    int numAtomsInCell = devCellListCount_[cellIndex];
                    for(size_t i=0; i<(size_t)numAtomsInCell; i++)
                    {
                        size_t cellListEntryIndex = CFunctorGenCellList::cellIndexToFlatIndex(Ix, Iy, Iz, i, maxAtomsInCell_, (size_t)cellCountX_, (size_t)cellCountY_);
                        if((int)cellListEntryIndex < maxCellListSize)
                        {
                            int atomIndex = devCellList_[cellListEntryIndex];
                            if(atomIndex < numAtoms_)
                            {
                                double drSqr = r0.distToAcrossPBC2(devAtomList_[atomIndex].r_, pbc_);

                                int currentNeigh = devNeighListCount_[atom.index_];
                                if((currentNeigh < maxNeighbors_) && (drSqr < (double)rCutoff2_))
                                {
                                    if(atomIndex == atom.index_) continue;
                                    devNeighList_[neighIndexToFlatIndex(atom.index_, currentNeigh, maxNeighbors_)] = atomIndex;
                                    devNeighListCount_[atom.index_]++;
                                }
                            }
                            else
                            {
                                return -4;
                            }
                        }
                        else
                        {
                            return -3;
                        }
                    }
                }
                else
                {
                    return -2;
                }
            }
        }
    }
    fflush(stdout);

    return devNeighListCount_[atom.index_];
}

HOSTDEV_CALLABLE size_t CFunctorGenNeighList::neighIndexToFlatIndex(size_t atomIndex, size_t neighIndex, int maxNeighbors)
{
    return neighIndex + maxNeighbors*atomIndex;
}

END_CUDA_COMPATIBLE()
