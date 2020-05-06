#include "MolecularTools.h"

double CMolecularTools::calcDistance2(const CAtom* fromAtom, const CAtom* toAtom, int frame, const C3DRect& pbc, bool distAcrossPBC) const
{
    if(distAcrossPBC) return fromAtom->getDistanceTo2UsingPBC(toAtom, frame, pbc);

    return fromAtom->getDistanceTo2(toAtom, frame);
}

double CMolecularTools::measureDist(const CAtom* atom1, const CAtom* atom2, int frame, const C3DRect* pbc)
{
    C3DVector r1 = atom1->r_[frame];
    C3DVector r2 = atom2->r_[frame];

    if(pbc) r1.moveToSameSideOfPBCAsThis(r2, *pbc);

    return (r2 - r1).norm();
}

double CMolecularTools::measureAngle(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, int frame, const C3DRect* pbc)
{
    C3DVector r1 = atom1->r_[frame];
    C3DVector r2 = atom2->r_[frame];
    C3DVector r3 = atom3->r_[frame];

    if(pbc)
    {
        r2.moveToSameSideOfPBCAsThis(r1, *pbc);
        r2.moveToSameSideOfPBCAsThis(r3, *pbc);
    }

    C3DVector v1 = r1 - r2;
    C3DVector v2 = r3 - r2;

    return v1.angle(v2);
}

double CMolecularTools::measureDihedral(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, const CAtom* atom4, int frame, const C3DRect* pbc)
{
    C3DVector r1 = atom1->r_[frame];
    C3DVector r2 = atom2->r_[frame];
    C3DVector r3 = atom3->r_[frame];
    C3DVector r4 = atom4->r_[frame];

    if(pbc)
    {
        r1.moveToSameSideOfPBCAsThis(r2, *pbc);
        r1.moveToSameSideOfPBCAsThis(r3, *pbc);
        r1.moveToSameSideOfPBCAsThis(r4, *pbc);
    }

    C3DVector a = r2 - r1;
    C3DVector b = r3 - r2;
    C3DVector c = r4 - r3;
    C3DVector n1 = a.cross(b);
    C3DVector n2 = b.cross(c);
    C3DVector k = n2.cross(b);
    double phi = n1.angle(n2);

    if((k * n1) < 0.0) phi = 2.0*M_PI - phi;

    return phi;
}

C3DVector CMolecularTools::rotatePosAroundBasisW(C3DVector basisPos, C3DBasis basis, C3DVector posToRotate, double relAngleToRotate)
{
    const double zeroWithin = 1.0E-6;
    C3DVector r_p = posToRotate - basisPos;

    // Project r_p into the new basis
    C3DVector R = C3DVector(r_p*basis.u_, r_p*basis.v_, r_p*basis.w_);
    C3DVector R_uv = R; R_uv.z_ = 0.0;
    double R_uv_abs = R_uv.norm();
    if(R_uv_abs <= zeroWithin) return posToRotate;

    // Rotate (Rx, Ry) around Basis.w
    double theta = atan2(R.y_, R.x_);

    C3DVector Rnew = C3DVector(R_uv_abs*cos(theta + relAngleToRotate),
                               R_uv_abs*sin(theta + relAngleToRotate),
                               R.z_);

    // Calculate r_new
    C3DVector r_new = basisPos + basis.u_*Rnew.x_ + basis.v_*Rnew.y_ + basis.w_*Rnew.z_;

    return r_new;
}

void CMolecularTools::modBondLengthTo(const CAtom* atom1, CAtom* atom2, double dist, int frame)
{
    int bondDest = atom2->getBondDestIndex(atom1);
    C3DVector newDist, oldDist, toMove;

    std::shared_ptr<std::vector<CAtom*>> atomsToMove = atom2->searchForLeafAtomsConnToBond(bondDest);

    newDist = atom2->r_[frame] - atom1->r_[frame];
    oldDist = newDist;
    newDist.normalize();
    newDist*= dist;

    toMove = newDist - oldDist;

    atom2->r_[frame]+= toMove;
    for(int i=0; i<atomsToMove->size(); i++)
    {
        (*atomsToMove)[i]->r_[frame]+= toMove;
    }
}

void CMolecularTools::modAngleTo(const CAtom* atom1, const CAtom* atom2, CAtom* atom3, double angle, int frame)
{
    int bondDest = atom3->getBondDestIndex(atom2);
    double dN1, R;
    C3DVector n1, n2, l, a, b;

    std::shared_ptr<std::vector<CAtom*>> atomsToMove = atom3->searchForLeafAtomsConnToBond(bondDest);

    a = atom2->r_[frame] - atom1->r_[frame];
    b = atom3->r_[frame] - atom2->r_[frame];

    // Find unit vectors n2, l, which are basis for
    // plane normal to b
    n1 = b.cross(a);
    dN1 = n1.norm();
    if(dN1 != 0.0)  n1*= (1.0 / dN1);
    else
    {
        // In this case, b and a are parallel and we need to find
        // a vector that is not parallel to a to generate n1. We
        // will always find that one of x1, x2, or x3 are non-
        // parallel!!!
        C3DVector x1(1.0, 0.0, 0.0);
        C3DVector x2(0.0, 1.0, 0.0);
        C3DVector x3(0.0, 0.0, 1.0);

        n1 = x1.cross(a);
        dN1 = n1.norm();
        if(dN1 != 0.0) n1*= (1.0 / dN1);
        else
        {
            n1 = x2.cross(a);
            dN1 = n1.norm();
            if(dN1 != 0.0) n1*= (1.0 / dN1);
            else
            {
                n1 = x3.cross(a);
                dN1 = n1.norm();
                n1*= (1.0 / dN1);
            }
        }
    }

    n2 = a.cross(n1);
    n2.normalize();

    l = a*(-1.0);
    l.normalize();

    // Find vector with correct angle in plane
    // defined by (n2, l) plane
    R = b.norm();
    atom3->r_[frame] = atom2->r_[frame] + (l*cos(angle) + n2*sin(angle))*R;

    // Rotate atoms connected to atom 3
    if((R > 0.0) && (l.norm2() > 0.0))
    {
        // In this usual case, where the length of a and b are non-zero,
        // the angle is actually changed. The same change must
        // occur on connected atoms

        // Find how much the atoms are to be rotated (i.e. DeltaPhi)
        // in the (n2, l) plane
        double lb = l*b;
        double lb_R = lb / R;
        if(lb_R > 1.0) lb_R = 1.0;
        if(lb_R < -1.0) lb_R = -1.0;
        double phi = acos(lb_R);
        double deltaPhi = angle - phi;

        // Rotate atoms i by DeltaPhi
        for(int i=0; i<atomsToMove->size(); i++)
        {
            if((*atomsToMove)[i] == atom1) continue;

            C3DVector Ci = (*atomsToMove)[i]->r_[frame] - atom2->r_[frame];

            // Find the norm of Ci's projection
            // onto the (l, n2) plane
            double n2Ci = n2*Ci;
            double n2Ci2 = n2Ci*n2Ci;
            double lCi = l*Ci;
            double lCi2 = lCi*lCi;
            double absProjCi = sqrt(n2Ci2 + lCi2);

            // Rotate the i'th atom only if proj(Ci) has a length.
            // If not, either Ci is parallel to the axis of rotation
            // and should not be rotated, or |Ci| = 0, in which case
            // Ci is invariant under any rotation
            if(absProjCi > 0.0)
            {
                // Find current angle of Ci in (n2, l) plane
                double lCi_dAbsProjCi = lCi / absProjCi;
                if(lCi_dAbsProjCi > 1.0) lCi_dAbsProjCi = 1.0;
                if(lCi_dAbsProjCi < -1.0) lCi_dAbsProjCi = -1.0;
                double phii = acos(lCi_dAbsProjCi);

                if(n2Ci < 0.0) phii = 2.0*M_PI - phii;

                // Find new angle of Ci in (n2, l) plane
                double phiiNew = phii + deltaPhi;

                // Calculate the new pos of atom 3 (i.e. proj along n1 + proj onto (n2, l) plane)
                (*atomsToMove)[i]->r_[frame] = atom2->r_[frame] + n1*(n1*Ci) + (l*cos(phiiNew) + n2*sin(phiiNew))*absProjCi;
            }
        }
    }
}

void CMolecularTools::modDihedralTo(const CAtom* atom1, const CAtom* atom2, const CAtom* atom3, CAtom* atom4, double angle, int frame)
{
    int bondDest = atom3->getBondDestIndex(atom2);
    double dN1, dN2;
    C3DVector n1, n2, l, k, a, b, c;

    std::shared_ptr<std::vector<CAtom*>> atomsToMove = atom3->searchForLeafAtomsConnToBond(bondDest);

    a = atom2->r_[frame] - atom1->r_[frame];
    b = atom3->r_[frame] - atom2->r_[frame];
    c = atom4->r_[frame] - atom3->r_[frame];

    // Find unit vector n1, which is normal to
    // plane formed by a and b
    n1 = a.cross(b);
    dN1 = n1.norm();

    if(dN1 != 0.0)  n1*= (1.0 / dN1);
    else
    {
        // In this case, b and a are parallel and we need to find
        // a vector that is not parallel to a to generate n1. We
        // will always find that one of x1, x2, or x3 are non-
        // parallel!!!
        C3DVector x1(1.0, 0.0, 0.0);
        C3DVector x2(0.0, 1.0, 0.0);
        C3DVector x3(0.0, 0.0, 1.0);

        n1 = x1.cross(a);
        dN1 = n1.norm();
        if(dN1 != 0.0) n1*= (1.0 / dN1);
        else
        {
            n1 = x2.cross(a);
            dN1 = n1.norm();
            if(dN1 != 0.0) n1*= (1.0 / dN1);
            else
            {
                n1 = x3.cross(a);
                dN1 = n1.norm();
                n1*= (1.0 / dN1);
            }
        }
    }

    // Find unit vector n2, which is normal to
    // plane formed by b and c
    n2 = b.cross(c);
    dN2 = n2.norm();
    if(dN2 != 0.0)  n2*= (1.0 / dN2);
    else
    {
        // In this case, b and c are parallel and the dihedral
        // is not defined properly!!!
        printf("Error: The atom2->atom3 vector is parallel to the atom3->atom4 vector: Dihedral not defined!\r\n");
        return;
    }

    // Define the orthonormal coordinate system (n2, l, k)
    // where l is parallel to b
    l = b;
    l.normalize();

    k = n2.cross(l);

    // Find current dihedral angle, dPhi, and the angle to
    // add to the current dihedral angle, dDeltaPhi
    double n1n2 = n1*n2;
    if(n1n2 > 1.0) n1n2 = 1.0;
    if(n1n2 < -1.0) n1n2 = -1.0;
    double phi = acos(n1n2);
    if(k*n1 < 0.0) phi = 2.0*M_PI - phi;
    double deltaPhi = angle - phi;

    // Rotate atoms i by DeltaPhi around l
    for(int i=0; i<atomsToMove->size(); i++)
    {
        if((*atomsToMove)[i] == atom1) continue;
        if((*atomsToMove)[i] == atom2) continue;
        if((*atomsToMove)[i] == atom3) continue;

        C3DVector Ci = (*atomsToMove)[i]->r_[frame] - atom3->r_[frame];

        // Find the norm of Ci's projection
        // onto the (k, n2) plane
        double n2Ci = n2*Ci;
        double n2Ci2 = n2Ci*n2Ci;
        double kCi = k*Ci;
        double kCi2 = kCi*kCi;
        double absProjCi = sqrt(n2Ci2 + kCi2);

        // Rotate the i'th atom only if proj(Ci) has a length.
        // If not, either Ci is parallel to the axis of rotation
        // and should not be rotated, or |Ci| = 0, in which case
        // Ci is invariant under any rotation
        if(absProjCi > 0.0)
        {
            // Find current angle of Ci in (n2, k) plane
            double kCi_dAbsProjCi = kCi / absProjCi;
            if(kCi_dAbsProjCi > 1.0) kCi_dAbsProjCi = 1.0;
            if(kCi_dAbsProjCi < -1.0) kCi_dAbsProjCi = -1.0;
            double phii = acos(kCi_dAbsProjCi);

            if(n2Ci < 0.0) phii = 2.0*M_PI - phii;

            // Find new angle of Ci in (n2, k) plane
            double phiiNew = phii + deltaPhi;

            // Calculate the new pos of atom 3 (i.e. proj along n1 + proj onto (n2, l) plane)
            (*atomsToMove)[i]->r_[frame] = atom3->r_[frame] + l*(l*Ci) + (k*cos(phiiNew) + n2*sin(phiiNew))*absProjCi;
        }
    }
}
