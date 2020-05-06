#include <iostream>
#include <math.h>
#include <vector>
#include "CmdAtom.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdAtom::getCmd()
{
    return "atom";
}

std::vector<std::string> CCmdAtom::getCmdLineKeywords()
{
    return { "atom", "at", "from", "dist", "dir", "bond", "angle", "cubecpy", "spherecpy", "random" };
}

std::vector<std::string> CCmdAtom::getCmdHelpLines()
{
    return {
                "atom <ID> at <x> <y> <z> [cubecpy <nx> <ny> <nz> <dx> <dy> <dz>, spherecpy <N> <R> [random]]",
                "atom <ID> from atom <n> dist <d> [cubecpy...]",
                "atom <ID> from atom <n> dir <dx> <dy> <dz> dist <d> [cubecpy...]",
                "atom <ID> from bond <n1> <n2> angle <angle> <dih> dist <d> [cubecpy...]",
                "atom <ID> from angle <n1> <n2> <n3> angle <angle> <dih> dist <d> [cubecpy...]"
           };
}

std::string CCmdAtom::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tIn general, ID could for example be C, H, or O (i.e. carbon, hydrogen, or\r\n";
    text+= "\toxygen). It is also possible to give names such as C1 or C2. As long as\r\n";
    text+= "\tthe name contains a C it is recognized as a carbon atom. Similarly, any\r\n";
    text+= "\tname containing O will be recognized as oxygen, etc. \r\n";
    text+= "\r\n";

    text+= "\tWhen atoms are added they attain a new index starting at index 0. The list\r\n";
    text+= "\tof atomic indices are obtained through the 'list' command.\r\n";

    text+= "\tVersions of the add command:\r\n";
    text+= "\r\n";
    text+= "\t 1) Add atom at the coordinate (<x>, <y>, <z>). Grid copy makes <na> copies separated\r\n";
    text+= "\t    by a distance <da>, where a is x, y, and z. Sphere-copy places <N> atoms at surface\r\n";
    text+= "\t    of a sphere of radius R, approx. equidistant., unless 'random' is specified.\r\n";
    text+= "\r\n";
    text+= "\t 2) Add atom a distance <d> in the y-direction from the atom with index <n>.\r\n";
    text+= "\r\n";
    text+= "\t 3) Add atom a distance <d> in the (<dx>, <dy>, <dz>)-direction from the atom\r\n";
    text+= "\t    with index <n>.\r\n";
    text+= "\r\n";
    text+= "\t 4) Add atom with index n3, such that the angle <n1>-<n2>-<n3> is <angle>\r\n";
    text+= "\t    degrees, with the <n2>-<n3> bond being rotated <dih> around the <n1>-<n2>\r\n";
    text+= "\t    bond (automatically using the x-, y-, or z- axis as reference). The bond\r\n";
    text+= "\t    length <n2>-<n3> is set to <d>.\r\n";
    text+= "\r\n";
    text+= "\t 5) Add atom with index n4, such that the angle <n2>-<n3>-<n4> is <angle>\r\n";
    text+= "\t    degrees, with the <n1>-<n2>-<n3>-<n4> dihedral set to <dih> degrees.\r\n";
    text+= "\t    The bond length <n3>-<n4> is set to <d>.";

    return text;
}

std::string CCmdAtom::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CAtom atom;

    text = CASCIIUtility::getArg(arguments, arg++);
    atom.setID(text.data());
    atom.sigma_ = state_->defaultAtProp_.getWDWRadius(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "at")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        atom.r_[0].x_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atom.r_[0].y_ = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atom.r_[0].z_ = atof(text.data());

        text = CASCIIUtility::getArg(arguments, arg++);
        if(text == "cubecpy")
        {
            addAtomByCubeCopy(atom, arguments, arg);
        }
        else if(text == "spherecpy")
        {
            if(!addAtomBySphereCopy(atom, arguments, arg)) return lastError_;
        }
        else
        {
            state_->addAtom(atom);
        }
    }
    else if(text == "from")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text == "atom")
        {
            int atomIndex;
            CAtom* refAtom;

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex = (int)atof(text.data());

            if(atomIndex < state_->atoms_.size())
            {
                refAtom = state_->atoms_[atomIndex].get();
                if(refAtom)
                {
                    text = CASCIIUtility::getArg(arguments, arg++);
                    if(text == "dist")
                    {
                        double dist;

                        text = CASCIIUtility::getArg(arguments, arg++);
                        dist = atof(text.data());

                        if(state_->currentFrame_ < refAtom->r_.size())
                        {
                            atom.r_[0] = refAtom->r_[state_->currentFrame_] + C3DVector(0.0, dist, 0.0);
                            if(text == "cubecpy")
                            {
                                addAtomByCubeCopy(atom, arguments, arg);
                            }
                            else if(text == "spherecpy")
                            {
                                if(!addAtomBySphereCopy(atom, arguments, arg)) return lastError_;
                            }
                            else
                            {
                                state_->addAtom(atom);
                            }
                        }
                    }
                    else if(text == "dir")
                    {
                        C3DVector   vDir;

                        text = CASCIIUtility::getArg(arguments, arg++);
                        vDir.x_ = atof(text.data());
                        text = CASCIIUtility::getArg(arguments, arg++);
                        vDir.y_ = atof(text.data());
                        text = CASCIIUtility::getArg(arguments, arg++);
                        vDir.z_ = atof(text.data());

                        vDir.normalize();

                        text = CASCIIUtility::getArg(arguments, arg++);
                        if(text == "dist")
                        {
                            double  dist;

                            text = CASCIIUtility::getArg(arguments, arg++);
                            dist = atof(text.data());
                            vDir*= dist;

                            if(state_->currentFrame_ < refAtom->r_.size())
                            {
                                atom.r_[0] = refAtom->r_[state_->currentFrame_] + vDir;
                                if(text == "cubecpy")
                                {
                                    addAtomByCubeCopy(atom, arguments, arg);
                                }
                                else if(text == "spherecpy")
                                {
                                    if(!addAtomBySphereCopy(atom, arguments, arg)) return lastError_;
                                }
                                else
                                {
                                    state_->addAtom(atom);
                                }
                            }
                        }
                        else
                        {
                            lastError_ = "Syntax Error: Ninth argument should be 'dist!";
                        }
                    }
                    else
                    {
                        lastError_ = "Syntax Error: Seventh argument should be 'dist' or 'dir!";
                    }
                }
                else
                {
                    lastError_ = "Could not find requested atom!";
                }
            }
            else
            {
                lastError_ = "Invalid atom index!";
            }
        }
        else if(text == "bond")
        {
            int atomIndex1, atomIndex2;
            CAtom* refAtom1;
            CAtom* refAtom2;

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex1 = (int)atof(text.data());

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex2 = (int)atof(text.data());

            if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size()))
            {
                refAtom1 = state_->atoms_[atomIndex1].get();
                refAtom2 = state_->atoms_[atomIndex2].get();
                if(refAtom1 && refAtom2)
                {
                    text = CASCIIUtility::getArg(arguments, arg++);
                    if(text =="angle")
                    {
                        double angleAroundBond, angleDirBond;

                        text = CASCIIUtility::getArg(arguments, arg++);
                        angleDirBond = atof(text.data());

                        text = CASCIIUtility::getArg(arguments, arg++);
                        angleAroundBond = atof(text.data());

                        text = CASCIIUtility::getArg(arguments, arg++);
                        if(text == "dist")
                        {
                            double len;
                            C3DVector vDir;

                            text = CASCIIUtility::getArg(arguments, arg++);
                            len = atof(text.data());

                            if((state_->currentFrame_ < refAtom1->r_.size()) && (state_->currentFrame_ < refAtom2->r_.size()))
                            {
                                vDir = calcDirVecFromBond(refAtom1, refAtom2, angleAroundBond, angleDirBond, len, state_->currentFrame_);

                                atom.r_[0] = refAtom2->r_[state_->currentFrame_] + vDir;
                                if(text == "cubecpy")
                                {
                                    addAtomByCubeCopy(atom, arguments, arg);
                                }
                                else if(text == "spherecpy")
                                {
                                    if(!addAtomBySphereCopy(atom, arguments, arg)) return lastError_;
                                }
                                else
                                {
                                    state_->addAtom(atom);
                                }
                            }
                        }
                        else
                        {
                            lastError_ = "Syntax Error: Tenth argument should be 'dist'!";
                        }
                    }
                    else
                    {
                        lastError_ = "Syntax Error: Seventh argument should be 'angle'!";
                    }
                }
                else
                {
                    lastError_ = "Could not find requested atom!";
                }
            }
            else
            {
                lastError_ = "Invalid atom index!";
            }
        }
        else if(text == "angle")
        {
            int atomIndex1, atomIndex2, atomIndex3;
            CAtom* refAtom1;
            CAtom* refAtom2;
            CAtom* refAtom3;

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex1 = (int)atof(text.data());

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex2 = (int)atof(text.data());

            text = CASCIIUtility::getArg(arguments, arg++);
            atomIndex3 = (int)atof(text.data());

            if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size())
               && (atomIndex3 < state_->atoms_.size()))
            {
                refAtom1 = state_->atoms_[atomIndex1].get();
                refAtom2 = state_->atoms_[atomIndex2].get();
                refAtom3 = state_->atoms_[atomIndex3].get();
                if(refAtom1 && refAtom2 && refAtom3)
                {
                    text = CASCIIUtility::getArg(arguments, arg++);
                    if(text == "angle")
                    {
                        double angleAroundBond, angleDirBond;

                        text = CASCIIUtility::getArg(arguments, arg++);
                        angleDirBond = atof(text.data());

                        text = CASCIIUtility::getArg(arguments, arg++);
                        angleAroundBond = atof(text.data());

                        text = CASCIIUtility::getArg(arguments, arg++);
                        if(text == "dist")
                        {
                            bool dihNotDef;
                            double len;
                            C3DVector vDir;

                            text = CASCIIUtility::getArg(arguments, arg++);
                            len = atof(text.data());

                            if((state_->currentFrame_ < refAtom1->r_.size()) && (state_->currentFrame_ < refAtom2->r_.size()) && (state_->currentFrame_ < refAtom3->r_.size()))
                            {
                                vDir = calcDirVecFromAngle(refAtom1, refAtom2, refAtom3, angleAroundBond, angleDirBond, len, dihNotDef, state_->currentFrame_);
                                if(dihNotDef)
                                {
                                    vDir = calcDirVecFromBond(refAtom2, refAtom3, angleAroundBond, angleDirBond, len, state_->currentFrame_);
                                    printf("Warning: Dihedral not defined. Used arbitrary choice of reference for dihedral angle!");
                                }

                                atom.r_[0] = refAtom3->r_[state_->currentFrame_] + vDir;
                                if(text == "cubecpy")
                                {
                                    addAtomByCubeCopy(atom, arguments, arg);
                                }
                                else if(text == "spherecpy")
                                {
                                    if(!addAtomBySphereCopy(atom, arguments, arg)) return lastError_;
                                }
                                else
                                {
                                    state_->addAtom(atom);
                                }
                            }
                        }
                        else
                        {
                            lastError_ = "Syntax Error: Eleventh argument should be 'dist'!";
                        }
                    }
                    else
                    {
                        lastError_ = "Syntax Error: Eighth argument should be 'angle'!";
                    }
                }
                else
                {
                    lastError_ = "Could not find requested atom!";
                }
            }
            else
            {
                lastError_ = "Invalid atom index!";
            }
        }
        else
        {
            lastError_ = "Syntax Error: Fifth argument should be a valid 'from' reference (i.e. from where?)!";
        }
    }
    else
    {
        lastError_ = "Syntax Error: Fourth argument should be how to add atom!";
    }

    return lastError_;
}

void CCmdAtom::addAtomByCubeCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg)
{
    std::string text;

    int     nX, nY, nZ;
    double  dX, dY, dZ;
    double  sX, sY, sZ;

    text = CASCIIUtility::getArg(arguments, arg++);
    nX = atoi(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    nY = atoi(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    nZ = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    dX = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    dY = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    dZ = atof(text.data());

    sX = atom.r_[0].x_;
    sY = atom.r_[0].y_;
    sZ = atom.r_[0].z_;

    for(int iX=0; iX<nX; iX++)
    {
        for(int iY=0; iY<nY; iY++)
        {
            for(int iZ=0; iZ<nZ; iZ++)
            {
                CAtom atomCpy = atom;
                atomCpy.r_[0].x_ = sX + double(iX)*dX;
                atomCpy.r_[0].y_ = sY + double(iY)*dY;
                atomCpy.r_[0].z_ = sZ + double(iZ)*dZ;
                state_->addAtom(atomCpy);
            }
        }
    }
}

bool CCmdAtom::addAtomBySphereCopy(const CAtom& atom, const std::vector<std::string>& arguments, size_t& arg)
{
    std::string text;

    int N, nRefinements=2000;
    bool random = false;
    double R, sX, sY, sZ;

    text = CASCIIUtility::getArg(arguments, arg++);
    N = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    R = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text =="random") random = true;

    if(N < 1)
    {
        lastError_ = std::string("Error: Not possible to have less than N=1 copies (N=") + std::to_string(N) + std::string(")!");
        return false;
    }

    sX = atom.r_[0].x_;
    sY = atom.r_[0].y_;
    sZ = atom.r_[0].z_;

    //////////////////////////////////////////////
    // Distribute N atoms, more or less uniformly
    // at the surface of a sphere of radius R
    //////////////////////////////////////////////
    std::vector<C3DVector> pointsOnSphere;
    std::vector<C3DVector> nextPointsOnSphere;


    // Step 1: Randomly select N points on the sphere
    srand((unsigned int)time(nullptr));
    for(int i=0; i<N; i++)
    {
        double theta = (double(rand()) / double(RAND_MAX)) * 2.0 * M_PI;
        double phi = (double(rand()) / double(RAND_MAX)) * M_PI;
        C3DVector v(R*cos(theta)*sin(phi), R*sin(theta)*sin(phi), R*cos(phi));

        pointsOnSphere.emplace_back(v);
    }


    // Step 2: Use molecular dynamics to spread particles evenly on the surface
    if(!random)
    {
        for(int t=0; t<nRefinements; t++)
        {
            // Calculate next positions for all particles and store them
            nextPointsOnSphere.clear();
            for(int i=0; i<pointsOnSphere.size(); i++)
            {
                // Calculate force F_i = \sum_j 1/(r_i - r_j)^2 on paricle i,
                // which in this case comes from Coulomb like terms
                C3DVector F_i;

                for(int j=0; j<pointsOnSphere.size(); j++)
                {
                    C3DVector r_ij, r_ij_u;
                    double r_ij2, f_ij;

                    if(j == i) continue;

                    r_ij = pointsOnSphere[i] - pointsOnSphere[j];
                    r_ij_u = r_ij.unit();
                    r_ij2 = r_ij.norm2();

                    if(r_ij2 != 0.0) f_ij = 1.0 / r_ij2;
                    else             f_ij = 1.0 / 0.0001;

                    F_i.x_+= r_ij_u.x_ * f_ij;
                    F_i.y_+= r_ij_u.y_ * f_ij;
                    F_i.z_+= r_ij_u.z_ * f_ij;
                }

                // Calculate unit vectors of plane on surface, u and v
                C3DVector r_i = pointsOnSphere[i];
                double absXY = sqrt(r_i.x_*r_i.x_ + r_i.y_*r_i.y_);
                double absXYZ = r_i.norm();
                double cosTheta = r_i.x_ / ( (absXY != 0.0) ? absXY : 0.0001 );
                double sinTheta = r_i.y_ / ( (absXY != 0.0) ? absXY : 0.0001 );
                double cosPhi = r_i.z_ / ( (absXYZ != 0.0) ? absXYZ : 0.0001 );
                double sinPhi = absXY / ( (absXYZ != 0.0) ? absXYZ : 0.0001 );
                C3DVector dR_dtheta(-R*sinTheta*sinPhi, R*cosTheta*sinPhi, 0.0);
                C3DVector dR_dphi(R*cosTheta*cosPhi, R*sinTheta*cosPhi, -R*sinPhi);
                double ABSdR_dtheta = dR_dtheta.norm();
                double ABSdR_dphi = dR_dphi.norm();
                C3DVector u(dR_dtheta.x_ / ((ABSdR_dtheta != 0.0) ? ABSdR_dtheta : 0.0001), dR_dtheta.y_ / ((ABSdR_dtheta != 0.0) ? ABSdR_dtheta : 0.0001), dR_dtheta.z_ / ((ABSdR_dtheta != 0.0) ? ABSdR_dtheta : 0.0001));
                C3DVector v(dR_dphi.x_ / ((ABSdR_dphi != 0.0) ? ABSdR_dphi : 0.0001), dR_dphi.y_ / ((ABSdR_dphi != 0.0) ? ABSdR_dphi : 0.0001), dR_dphi.z_ / ((ABSdR_dphi != 0.0) ? ABSdR_dphi : 0.0001));
                double U = F_i * u;
                double V = F_i * v;
                C3DVector projF_i = u*U + v*V;

                // Perform molecular dynamics step
                projF_i.normalize();
                double epsilon = 2.0*M_PI*R / 1000.0;
                C3DVector r_i_tp1 = r_i + projF_i*epsilon;

                nextPointsOnSphere.emplace_back(r_i_tp1);
            }

            // Copy all next positions to current array of psitions
            for(int i=0; i<nextPointsOnSphere.size(); i++)
            {
                if(i < pointsOnSphere.size()) pointsOnSphere[i] = nextPointsOnSphere[i];
            }
        }
    }


    // Step 3: Place points on sphere at pos (sX,sY,sZ)
    for(int i=0; i<N; i++)
    {
        CAtom atomCpy = atom;
        atomCpy.r_[0].x_ = sX + pointsOnSphere[i].x_;
        atomCpy.r_[0].y_ = sY + pointsOnSphere[i].y_;
        atomCpy.r_[0].z_ = sZ + pointsOnSphere[i].z_;
        state_->addAtom(atomCpy);
    }

    return true;
}

C3DVector CCmdAtom::calcDirVecFromBond(const CAtom* at1, const CAtom* at2, double angleAroundBond, double angleDirBond, double len, int frame) const
{
    C3DVector dir;
    C3DVector nyz(1.0, 0.0, 0.0);
    C3DVector nzx(0.0, 1.0, 0.0);
    C3DVector nxy(0.0, 0.0, 1.0);
    C3DVector x1(0.0, 1.0, 0.0);
    C3DVector x2(0.0, 0.0, 1.0);
    C3DVector w, u, v, k, l;
    C3DVector bond = at2->r_[frame] - at1->r_[frame];
    double norm1, norm2;

    angleAroundBond*= M_PI / 180.0;
    angleDirBond*= M_PI / 180.0;

    // Find which planes yz, zx, or xy the bond vector is 'more' normal onto
    norm1 = (nyz.cross(bond)).norm2();
    norm2 = (nzx.cross(bond)).norm2();
    if(norm2 < norm1)
    {
        norm1 = norm2;
        x1 = C3DVector(0.0, 0.0, 1.0);
        x2 = C3DVector(1.0, 0.0, 0.0);
    }
    norm2 = (nxy.cross(bond)).norm2();
    if(norm2 < norm1)
    {
        norm1 = norm2;
        x1 = C3DVector(1.0, 0.0, 0.0);
        x2 = C3DVector(0.0, 1.0, 0.0);
    }

    // Find which plane basis x1 or x2 is 'less' parallel to bond vector
    norm1 = x1*bond;
    norm2 = x2*bond;
    if(norm2 < norm1) w = x2;
    else                w = x1;

    // Calculate basis u, v, for plane normal to bond vector
    u = bond.cross(w);
    u.normalize();

    v = u.cross(bond);
    v.normalize();

    // Calculate basis k, l, for plane containing bond vector
    k = u*cos(angleAroundBond) + v*sin(angleAroundBond);
    l = bond*(-1.0);
    l.normalize();

    // Calculate direction vector
    dir = l*cos(angleDirBond) + k*sin(angleDirBond);
    dir*= len;

    return dir;
}

C3DVector CCmdAtom::calcDirVecFromAngle(const CAtom* at1, const CAtom* at2, const CAtom* at3,
                                    double angleAroundBond, double angleDirBond, double len, bool& dihedralNotDefined, int frame) const
{
    C3DVector dir;
    C3DVector a, b, k, l;
    C3DVector n1, n2;
    double dn1, dn2;

    angleAroundBond*= M_PI / 180.0;
    angleDirBond*= M_PI / 180.0;

    dihedralNotDefined = false;

    // Define bond vectors between reference atoms
    a = at1->r_[frame] - at2->r_[frame];
    b = at3->r_[frame] - at2->r_[frame];

    // Calculate plane normal to b, defined by n1, and n2, where
    // angle is defined with respect to a
    n1 = b.cross(a);
    n2 = n1.cross(b);

    dn1 = n1.norm();
    dn2 = n2.norm();

    if((dn1 == 0.0) || (dn2 == 0.0))
    {
        dihedralNotDefined = true;
        return dir;
    }

    n1*= (1.0 / dn1);
    n2*= (1.0 / dn2);

    // Calculate unit vector k in (n1,n2)-plane, which defines
    // the angle around the b vector (i.e. dihedral angle)
    k = n2*cos(angleAroundBond) + n1*sin(angleAroundBond);

    // Calc vector l along b, that together with k forms a plane
    // with b in it, as well as the new bond to be created
    l = b*(-1.0);
    l.normalize();

    // Calculate the displacement of the new atom from the atom
    // at b
    dir = l*cos(angleDirBond) + k*sin(angleDirBond);
    dir*= len;

    return dir;
}
