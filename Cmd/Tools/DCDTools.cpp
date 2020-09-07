#include "DCDTools.h"

double CDCDTools::CHBondCriteria::calcSizeOfGridToFitHBondPairs() const
{
    std::vector<double> distances;

    // To estimate, we will find the longest lenght (not equal to DBL_MAX)
    // and double that length to estimate the minumum length of a d-h...a pair.
    // If all lenghts are DBL_MAX, then return a generic length of 3AA
    if(fabs((double)maxLenR_dh_ - DBL_MAX) > double(FLT_MIN) ) distances.emplace_back(maxLenR_dh_);
    if(fabs((double)maxLenR_ha_ - DBL_MAX) > double(FLT_MIN)) distances.emplace_back(maxLenR_ha_);
    if(fabs((double)maxLenR_da_ - DBL_MAX) > double(FLT_MIN)) distances.emplace_back(maxLenR_da_);

    double maxLen = 0.0;
    for(int i=0; i<distances.size(); i++)
    {
        if(distances[i] > maxLen) maxLen = distances[i];
    }

    return  (maxLen > 0.0) ? (maxLen + 0.1) : 3.0;
}

C3DVector CDCDTools::getMoleculeDipoleMoment(const std::vector<int>& atomIndices, const CDCDFile* dcdFile, C3DVector& Rc, bool chargeNeutralFormulation) const
{
    C3DVector P, Rc_int;

    Rc = getCenterOfMass(atomIndices, dcdFile);
    if(!chargeNeutralFormulation) Rc_int = Rc;

    for(int i=0; i<atomIndices.size(); i++)
    {
        double q_i = state_->atoms_[atomIndices[i]]->Q_;
        C3DVector r_i = dcdFile->getCoordinate(atomIndices[i]);
        C3DVector d_i = r_i - Rc_int;
        P+= (d_i * q_i);
    }

    return P;
}

void CDCDTools::atomicUnwrapCurrDCDFrame(CDCDFile* dcdFile, const std::vector<std::vector<int>>& molList, const C3DRect* pbc) const
{
    if(!dcdFile) return;

    C3DRect PBC;

    if(pbc) PBC = *pbc;
    else PBC = dcdFile->getCurrentPBC();

    C3DVector pbcWidth_2(PBC.getWidthX() / 2.0, PBC.getWidthY() / 2.0, PBC.getWidthZ() / 2.0);
    C3DVector pbcWidth(PBC.getWidthX(), PBC.getWidthY(), PBC.getWidthZ());
    int numAtomsInRec = dcdFile->getNumCoordinatesInRecord();
    for(int i=0; i<molList.size(); i++)
    {
        // First check if the molecule needs to be unwrapped in any directions
        bool needUnwrap[3] = { false, false, false };
        for(int j=0; j<molList[i].size(); j++)
        {
            if(molList[i][j] >= numAtomsInRec) continue;
            C3DVector r_j = dcdFile->getCoordinate(molList[i][j]);

            for(int k=j+1; k<molList[i].size(); k++)
            {
                if(molList[i][k] >= numAtomsInRec) continue;
                C3DVector r_k = dcdFile->getCoordinate(molList[i][k]);

                C3DVector R_jk = r_j - r_k;
                if(fabs(R_jk.x_) > pbcWidth_2.x_) needUnwrap[0] = true;
                if(fabs(R_jk.y_) > pbcWidth_2.y_) needUnwrap[1] = true;
                if(fabs(R_jk.z_) > pbcWidth_2.z_) needUnwrap[2] = true;
            }
        }

        // Unwrap a single molecule. Always move atoms located in the lower part
        // of PBCs (below PBC/2, in x,y, or z direction) to higher part of PBCs
        if(!needUnwrap[0] && !needUnwrap[1] && !needUnwrap[2]) continue;
        for(int j=0; j<molList[i].size(); j++)
        {
            if(molList[i][j] >= numAtomsInRec) continue;
            C3DVector r = dcdFile->getCoordinate(molList[i][j]);
            C3DVector dist = PBC.rHigh_ - r;

            // If unwrap is needed in x-direction
            if(needUnwrap[0])
            {
                if(dist.x_ > pbcWidth_2.x_) dcdFile->setCoordinate(molList[i][j], 0, r.x_ + pbcWidth.x_);
            }
            // If unwrap is needed in y-direction
            if(needUnwrap[1])
            {
                if(dist.y_ > pbcWidth_2.y_) dcdFile->setCoordinate(molList[i][j], 1, r.y_ + pbcWidth.y_);
            }
            // If unwrap is needed in z-direction
            if(needUnwrap[2])
            {
                if(dist.z_ > pbcWidth_2.z_) dcdFile->setCoordinate(molList[i][j], 2, r.z_ + pbcWidth.z_);
            }
        }
    }
}

C3DVector CDCDTools::getCenterOfMass(const std::vector<int>& atomIndices, const CDCDFile* dcdFile) const
{
    double M = 0.0;
    C3DVector Rc;

    for(int i=0; i<atomIndices.size(); i++)
    {
        double m = state_->atoms_[atomIndices[i]]->m_;
        if(fabs(m) < 1E-5)
        {
            printf("Error: found zero mass!");
            return Rc;
        }
        M+= m;

        Rc+= dcdFile->getCoordinate(atomIndices[i])*m;
    }
    if(M != 0.0) Rc*= (1.0 / M);
    else
    {
        printf("Error: found vanishing total mass!");
    }

    return Rc;
}

std::shared_ptr<std::vector<CDCDTools::CHBond>> CDCDTools::retrieveHBondsCurrDCDFrame(const CDCDFile* dcdFile, const C3DRect* pbc, const std::vector<CHBondCriteria>& HBondCriteria, EHBondSpan HBondSpan, bool noPBC) const
{
    int numAtoms = dcdFile->getNumCoordinatesInRecord();
    C3DRect pbcOuter, pbcRim;


    if(pbc) pbcOuter = *pbc;
    else pbcOuter = dcdFile->getCurrentPBC();

    pbcRim = dcdFile->getCurrentBoundingBox();
    pbcRim.expandByLength(0.1);

    auto HBonds = std::make_shared<std::vector<CHBond>>();
    for(int i=0; i<HBondCriteria.size(); i++)
    {
        // Calculate grid size and number of grid points
        double gridSize = HBondCriteria[i].calcSizeOfGridToFitHBondPairs();
        int N[3];
        N[0] = (int)ceil(pbcRim.getWidthX() / gridSize);
        N[1] = (int)ceil(pbcRim.getWidthY() / gridSize);
        N[2] = (int)ceil(pbcRim.getWidthZ() / gridSize);


        // Set up donor, hydrogen and acceptor grids. Each grid point contains an
        // array of atom indices.
        std::vector<std::vector<std::vector<std::vector<int>>>> a3DGridOfDonors;
        a3DGridOfDonors.resize(N[0]);
        for(int j=0; j<a3DGridOfDonors.size(); j++)
        {
            a3DGridOfDonors[j].resize(N[1]);
            for(int k=0; k<a3DGridOfDonors[j].size(); k++)
                a3DGridOfDonors[j][k].resize(N[2]);
        }

        std::vector<std::vector<std::vector<std::vector<int>>>> a3DGridOfHydrogens;
        a3DGridOfHydrogens.resize(N[0]);
        for(int j=0; j<a3DGridOfHydrogens.size(); j++)
        {
            a3DGridOfHydrogens[j].resize(N[1]);
            for(int k=0; k<a3DGridOfHydrogens[j].size(); k++)
                a3DGridOfHydrogens[j][k].resize(N[2]);
        }

        std::vector<std::vector<std::vector<std::vector<int>>>> a3DGridOfAcceptors;
        a3DGridOfAcceptors.resize(N[0]);
        for(int j=0; j<a3DGridOfAcceptors.size(); j++)
        {
            a3DGridOfAcceptors[j].resize(N[1]);
            for(int k=0; k<a3DGridOfAcceptors[j].size(); k++)
                a3DGridOfAcceptors[j][k].resize(N[2]);
        }


        // Retrieve donor, hydrogen and acceptor indices, based on atom list data
        std::vector<int> donors, hydrogens, acceptors;
        for(int j=0; j<numAtoms; j++)
        {
            if(j >= state_->atoms_.size()) continue;
            std::string ID = state_->atoms_[j]->getID();
            if(ID == HBondCriteria[i].getDonor()) donors.emplace_back(j);
            if(ID == HBondCriteria[i].getHydrogen()) hydrogens.emplace_back(j);
            if(ID == HBondCriteria[i].getAcceptor()) acceptors.emplace_back(j);
        }


        // Bin donor, hydrogen and acceptor indices
        double dWidthRim[3] = { pbcRim.getWidthX(), pbcRim.getWidthY(), pbcRim.getWidthZ() };

        for(int j=0; j<donors.size(); j++)
        {
            C3DVector r = dcdFile->getCoordinate(donors[j]);
            int iBinIndex[3];
            for(int k=0; k<3; k++)
                iBinIndex[k] = (int)((r[k] - pbcRim.rLow_[k]) / dWidthRim[k] * double(N[k]));

            if((iBinIndex[0] < 0) || (iBinIndex[0] >= N[0])) continue;
            if((iBinIndex[1] < 0) || (iBinIndex[1] >= N[1])) continue;
            if((iBinIndex[2] < 0) || (iBinIndex[2] >= N[2])) continue;

            a3DGridOfDonors[iBinIndex[0]][iBinIndex[1]][iBinIndex[2]].emplace_back(donors[j]);
        }

        for(int j=0; j<hydrogens.size(); j++)
        {
            C3DVector r = dcdFile->getCoordinate(hydrogens[j]);
            int iBinIndex[3];
            for(int k=0; k<3; k++)
                iBinIndex[k] = (int)((r[k] - pbcRim.rLow_[k]) / dWidthRim[k] * double(N[k]));

            if((iBinIndex[0] < 0) || (iBinIndex[0] >= N[0])) continue;
            if((iBinIndex[1] < 0) || (iBinIndex[1] >= N[1])) continue;
            if((iBinIndex[2] < 0) || (iBinIndex[2] >= N[2])) continue;

            a3DGridOfHydrogens[iBinIndex[0]][iBinIndex[1]][iBinIndex[2]].emplace_back(hydrogens[j]);
        }

        for(int j=0; j<acceptors.size(); j++)
        {
            C3DVector r = dcdFile->getCoordinate(acceptors[j]);
            int iBinIndex[3];
            for(int k=0; k<3; k++)
                iBinIndex[k] = (int)((r[k] - pbcRim.rLow_[k]) / dWidthRim[k] * double(N[k]));

            if((iBinIndex[0] < 0) || (iBinIndex[0] >= N[0])) continue;
            if((iBinIndex[1] < 0) || (iBinIndex[1] >= N[1])) continue;
            if((iBinIndex[2] < 0) || (iBinIndex[2] >= N[2])) continue;

            a3DGridOfAcceptors[iBinIndex[0]][iBinIndex[1]][iBinIndex[2]].emplace_back(acceptors[j]);
        }


        // Find list of potential hydrogen bonds (based on bond and h-bond lenght criteria alone)
        double R_dh_c2 = HBondCriteria[i].R_dh() * HBondCriteria[i].R_dh();
        double R_ha_c2 = HBondCriteria[i].R_ha() * HBondCriteria[i].R_ha();
        double R_da_c2 = HBondCriteria[i].R_da() * HBondCriteria[i].R_da();
        for(int j1=0; j1<a3DGridOfDonors.size(); j1++)
        {
            for(int j2=0; j2<a3DGridOfDonors[j1].size(); j2++)
            {
                for(int j3=0; j3<a3DGridOfDonors[j1][j2].size(); j3++)
                {
                    // Go through each donor atom within cell (j1, j2, j3)
                    for(int j4=0; j4<a3DGridOfDonors[j1][j2][j3].size(); j4++)
                    {
                        int iDonor = a3DGridOfDonors[j1][j2][j3][j4];
                        C3DVector rD = dcdFile->getCoordinate(iDonor);

                        // For each donor cell, go through each neighbor cell to search
                        // for hydrogen atoms, as well as in the same cell
                        for(int k1=(j1-1); k1<=(j1+1); k1++)
                        {
                            int k1pbc = k1;
                            if(!noPBC && (k1 < 0)) k1pbc+= N[0];
                            if(!noPBC && (k1 >= N[0])) k1pbc-= N[0];
                            if((k1pbc < 0) || (k1pbc >= N[0])) continue;

                            for(int k2=(j2-1); k2<=(j2+1); k2++)
                            {
                                int k2pbc = k2;
                                if(!noPBC && (k2 < 0)) k2pbc+= N[1];
                                if(!noPBC && (k2 >= N[1])) k2pbc-= N[1];
                                if((k2pbc < 0) || (k2pbc >= N[1])) continue;

                                for(int k3=(j3-1); k3<=(j3+1); k3++)
                                {
                                    int k3pbc = k3;
                                    if(!noPBC && (k3 < 0)) k3pbc+= N[2];
                                    if(!noPBC && (k3 >= N[2])) k3pbc-= N[2];
                                    if((k3pbc < 0) || (k3pbc >= N[2])) continue;

                                    // Go through each hydrogen atom within cell (k1, k2, k3)
                                    for(int k4=0; k4<a3DGridOfHydrogens[k1pbc][k2pbc][k3pbc].size(); k4++)
                                    {
                                        int indexHydrogen = a3DGridOfHydrogens[k1pbc][k2pbc][k3pbc][k4];
                                        C3DVector rH = dcdFile->getCoordinate(indexHydrogen);
                                        double R_dh2 = rD.distToAcrossPBC2(rH, pbcOuter);

                                        // Check if the bond criteria r_dh < r_dh^c is satisfied between
                                        // atom j4 in cell (j1, j2, j3) and atom k4 in this cell
                                        if(R_dh2 < R_dh_c2)
                                        {
                                            // For each detected hydrogen, go through each neighbor cell to search
                                            // for acceptor atoms, as well as in the same cell
                                            for(int l1=(k1pbc-1); l1<=(k1pbc+1); l1++)
                                            {
                                                int l1pbc = l1;
                                                if(!noPBC && (l1 < 0)) l1pbc+= N[0];
                                                if(!noPBC && (l1 >= N[0])) l1pbc-= N[0];
                                                if((l1pbc < 0) || (l1pbc >= N[0])) continue;

                                                for(int l2=(k2pbc-1); l2<=(k2pbc+1); l2++)
                                                {
                                                    int l2pbc = l2;
                                                    if(!noPBC && (l2 < 0)) l2pbc+= N[1];
                                                    if(!noPBC && (l2 >= N[1])) l2pbc-= N[1];
                                                    if((l2pbc < 0) || (l2pbc >= N[1])) continue;

                                                    for(int l3=(k3pbc-1); l3<=(k3pbc+1); l3++)
                                                    {
                                                        int l3pbc = l3;
                                                        if(!noPBC && (l3 < 0)) l3pbc+= N[2];
                                                        if(!noPBC && (l3 >= N[2])) l3pbc-= N[2];
                                                        if((l3pbc < 0) || (l3pbc >= N[2])) continue;

                                                        // Go through each acceptor atom within cell (l1, l2, l3)
                                                        for(int l4=0; l4<a3DGridOfAcceptors[l1pbc][l2pbc][l3pbc].size(); l4++)
                                                        {
                                                            int indexAcceptor = a3DGridOfAcceptors[l1pbc][l2pbc][l3pbc][l4];
                                                            C3DVector rA = dcdFile->getCoordinate(indexAcceptor);
                                                            double R_ha2 = rH.distToAcrossPBC2(rA, pbcOuter);

                                                            // Check if the bond criteria r_ha < r_ha^c is satisfied between
                                                            // atom k4 in cell (k1, k2, k3) and atom l4 in this cell
                                                            if((R_ha2 < R_ha_c2) && (iDonor != indexAcceptor))
                                                            {
                                                                // We found a d-h...a H-bond set satisfying
                                                                // r_dh < r_dh^c and r_ha < r_ha^c. Store!
                                                                CHBond HBond(iDonor, indexHydrogen, indexAcceptor);
                                                                HBonds->emplace_back(HBond);
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        // Remove entries in H-bond list not sattisfying H-bond criteria.
        // We have already tested r_dh and r_ha. This leaves r_da and
        // the d-h...a angle left for testing
        for(int j=0; j<HBonds->size(); j++)
        {
            C3DVector rD = dcdFile->getCoordinate((*HBonds)[j].d_);
            C3DVector rA = dcdFile->getCoordinate((*HBonds)[j].a_);
            double R_da2 = rD.distToAcrossPBC2(rA, pbcOuter);

            if(R_da2 >= R_da_c2)
            {
                HBonds->erase(HBonds->begin() + j);
                j--;
            }
        }

        double Alpha_dha_c = HBondCriteria[i].alpha_dha() * M_PI / 180.0;
        for(int j=0; j<HBonds->size(); j++)
        {
            C3DVector rD = dcdFile->getCoordinate((*HBonds)[j].d_);
            C3DVector rH = dcdFile->getCoordinate((*HBonds)[j].h_);
            C3DVector rA = dcdFile->getCoordinate((*HBonds)[j].a_);
            double alpha_dha = rH.angleAcrossPBC(rD, rA, pbcOuter);

            if(alpha_dha <= Alpha_dha_c)
            {
                HBonds->erase(HBonds->begin() + j);
                j--;
            }
        }


        // Remove intra-molecular H-bonds if neccessary. An H-bond is
        // considered intra-molecular if the atoms involved have the same
        // molecule ID within the list of atoms shown graphically in
        // the 3D View
        if(HBondSpan == spanOnlyInterMolecularHBonds)
        {
            for(int j=0; j<HBonds->size(); j++)
            {
                if((*HBonds)[j].d_ >= state_->atoms_.size()) continue;
                if((*HBonds)[j].a_ >= state_->atoms_.size()) continue;

                int donorMolIndex = state_->atoms_[(*HBonds)[j].d_]->getMolIndex();
                int acceptorMolIndex = state_->atoms_[(*HBonds)[j].a_]->getMolIndex();

                if(donorMolIndex == acceptorMolIndex)
                {
                    HBonds->erase(HBonds->begin() + j);
                    j--;
                }
            }
        }
        if(HBondSpan == spanOnlyIntraMolecularHBonds)
        {
            for(int j=0; j<HBonds->size(); j++)
            {
                if((*HBonds)[j].d_ >= state_->atoms_.size()) continue;
                if((*HBonds)[j].a_ >= state_->atoms_.size()) continue;

                int donorMolIndex = state_->atoms_[(*HBonds)[j].d_]->getMolIndex();
                int acceptorMolIndex = state_->atoms_[(*HBonds)[j].a_]->getMolIndex();

                if(donorMolIndex != acceptorMolIndex)
                {
                    HBonds->erase(HBonds->begin() + j);
                    j--;
                }
            }
        }
    }

    return HBonds;
}

double CDCDTools::measureDist(const CDCDFile& dcdFile, int index1, int index2, const C3DRect* pbc)
{
    C3DVector r1 = dcdFile.getCoordinate(index1);
    C3DVector r2 = dcdFile.getCoordinate(index2);

    if(pbc) r1.moveToSameSideOfPBCAsThis(r2, *pbc);

    return (r2 - r1).norm();
}

double CDCDTools::measureAngle(const CDCDFile& dcdFile, int index1, int index2, int index3, const C3DRect* pbc)
{
    C3DVector r1 = dcdFile.getCoordinate(index1);
    C3DVector r2 = dcdFile.getCoordinate(index2);
    C3DVector r3 = dcdFile.getCoordinate(index3);

    if(pbc)
    {
        r2.moveToSameSideOfPBCAsThis(r1, *pbc);
        r2.moveToSameSideOfPBCAsThis(r3, *pbc);
    }

    C3DVector v1 = r1 - r2;
    C3DVector v2 = r3 - r2;

    return v1.angle(v2);
}

double CDCDTools::measureDihedral(const CDCDFile& dcdFile, int index1, int index2, int index3, int index4, const C3DRect* pbc)
{
    C3DVector r1 = dcdFile.getCoordinate(index1);
    C3DVector r2 = dcdFile.getCoordinate(index2);
    C3DVector r3 = dcdFile.getCoordinate(index3);
    C3DVector r4 = dcdFile.getCoordinate(index4);

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
