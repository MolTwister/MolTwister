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
}

void CFunctorGenNeighList::setForceFieldMatrices(CMDFFMatrices& ffMatrices)
{
    devCellList_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellList_[0]);
    devCellListCount_ = mtraw_pointer_cast(&ffMatrices.cellList_.devCellListCount_[0]);
    devAtomCellIndicesRaw_ = mtraw_pointer_cast(&ffMatrices.cellList_.devAtomCellIndices_[0]);

    devNeighList_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighList_[0]);
    devNeighListCount_ = mtraw_pointer_cast(&ffMatrices.neighList_.devNeighListCount_[0]);

    cellCountX_ = ffMatrices.cellList_.getCellCountX();
    cellCountY_ = ffMatrices.cellList_.getCellCountY();
    cellCountZ_ = ffMatrices.cellList_.getCellCountZ();

    maxNeighbors_ = ffMatrices.neighList_.getMaxNeighbors();
    maxAtomsInCell_ = ffMatrices.cellList_.getMaxAtomsInCell();
}

HOSTDEV_CALLABLE int CFunctorGenNeighList::operator()(CMDFFMatrices::CAtom& atom)
{
    CMDFFMatrices::CCellListIndex cellIndex = devAtomCellIndicesRaw_[atom.index_];

    int ixl = cellIndex.ix_ - 1;
    int ixh = ixl + 2;

    int iyl = cellIndex.iy_ - 1;
    int iyh = iyl + 2;

    int izl = cellIndex.iz_ - 1;
    int izh = izl + 2;

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
                int numAtomsInCell = devCellListCount_[CFunctorGenCellList::cellIndexToFlatIndex(Ix, Iy, Iz, (size_t)cellCountX_, (size_t)cellCountY_)];
                for(size_t i=0; i<(size_t)numAtomsInCell; i++)
                {
                    int atomIndex = devCellList_[CFunctorGenCellList::cellIndexToFlatIndex(Ix, Iy, Iz, i, (size_t)maxAtomsInCell_, (size_t)cellCountX_, (size_t)cellCountY_)];
                    int currentNeigh = devNeighListCount_[atom.index_];
                    if(currentNeigh < maxNeighbors_)
                    {
                        devNeighList_[neighIndexToFlatIndex(atom.index_, currentNeigh, maxNeighbors_)] = atomIndex;
                        devNeighListCount_[atom.index_]++;
                    }
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

END_CUDA_COMPATIBLE()
