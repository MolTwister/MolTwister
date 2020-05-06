#pragma once
#include "ToolsBase.h"
#include "../../MolTwisterAtom.h"
#include "../../Utilities/3DRect.h"
#include "../../Utilities/3DBasis.h"

class CMolecularTools : public CToolsBase
{
public:
    CMolecularTools() = delete;
    CMolecularTools(CMolTwisterState* state, FILE* stdOut) : CToolsBase(state, stdOut) { }

public:
    static double measureDihedral(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, const CAtom* atom4, int frame, const C3DRect* pbc=nullptr);
    static double measureAngle(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, int frame, const C3DRect* pbc=nullptr);
    static double measureDist(const CAtom* atom1, const CAtom* atom2, int frame, const C3DRect* pbc=nullptr);
    double calcDistance2(const CAtom* fromAtom, const CAtom* toAtom, int frame, const C3DRect& pbc, bool distAcrossPBC) const;
    static C3DVector rotatePosAroundBasisW(C3DVector basisPos, C3DBasis basis, C3DVector posToRotate, double relAngleToRotate);
    static void modBondLengthTo(const CAtom* atom1, CAtom* atom2, double dist, int frame);
    static void modAngleTo(const CAtom* atom1, const CAtom* atom2, CAtom* atom3, double angle, int frame);
    static void modDihedralTo(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, CAtom* atom4, double angle, int frame);
};
