//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include <iostream>
#include <vector>
#include <ctime>
#include "Tools/MolTwisterStateTools.h"
#include "Tools/ProgressBar.h"
#include "MolTwisterCmdCopy.h"

void CCmdCopy::onAddKeywords()
{
    addKeyword("copy");
    addKeyword("atoms");
    addKeyword("names");
    addKeyword("resnames");
    addKeyword("all");
    addKeyword("sel");
    addKeyword("random");
}

std::string CCmdCopy::getHelpString() const
{
    std::string text;
    
    text+= "\tUsage: copy <specifier>\r\n";
    text+= "\r\n";
    text+= "\tMake several copies of atoms, or collection of atoms. The allowed specifiers are:\r\n";
    text+= "\r\n";
    text+= "\t        * atoms <nx> <ny> <nz> <dx> <dy> <dz> <index1> ... <indexN>         : Copy <ni> times\r\n";
    text+= "\t          with distance <di> between each copy, where i is x,y or z.\r\n";
    text+= "\t        * names <nx> <ny> <nz> <dx> <dy> <dz> <name1> ... <nameN>           : Copy based on names.\r\n";
    text+= "\t        * names random <len x> <len y> <len z> <n> <dr> <name1> ... <nameN> : Copy based on names.\r\n";
    text+= "\t        * all <nx> <ny> <nz> <dx> <dy> <dz>                                 : Copy all.\r\n";
    text+= "\t        * all random <len x> <len y> <len z> <n> <dr>                       : Copy all.\r\n";
    text+= "\t        * sel <nx> <ny> <nz> <dx> <dy> <dz>                                 : Copy only selected.\r\n";
    text+= "\t        * sel random <len x> <len y> <len z> <n> <dr>                       : Copy only selected.\r\n";
    text+= "\t        * resnames <nx> <ny> <nz> <dx> <dy> <dz> <resname1> ... <resnameN>  : Copy based on resnames.\r\n";
    text+= "\t        * resnames random <len x> <len y> <len z> <n> <dr> <resname1> ... <resnameN> : Copy based on resnames.\r\n";
    text+= "\r\n";
    text+= "\tIf the random keyword is used, then <n> copies of the atomic set are randomly placed and oriented within a box\r\n";
    text+= "\twith size <len x>, <len y>, and <len z>, starting from the first copy, such that the atoms of the sets\r\n";
    text+= "\tnever touch, with a least distance <dr>.";

    return text;
}

void CCmdCopy::execute(std::string commandLine)
{
    std::vector<std::string> objectsToCopy;
    std::vector<CAtom*> atomsToCopy;
    std::string text, commandString;
    double dX = 0.0, dY = 0.0, dZ = 0.0;
    int arg = 1, length, count = 0;
    int nX = 0, nY = 0, nZ = 0;
    double sX = 0.0, sY = 0.0, sZ = 0.0;
    double lenX = 0.0, lenY = 0.0, lenZ = 0.0;
    int N = 0;
    double dr = 0.0;
    
    if(!state_) return;
    
    commandString = CASCIIUtility::getWord(commandLine, arg++);
    if((commandString != "atoms") &&
       (commandString != "names") &&
       (commandString != "all") &&
       (commandString != "sel") &&
       (commandString != "resnames"))
    {
        printf("Syntax Error: First argument should be the type of object to copy!");
        return;
    }
    
    // Get amoount of copies to do and distances
    bool useRandomDistribution = false;
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text != "random")
    {
        nX = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        nY = atoi(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        nZ = atoi(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        dX = atof(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        dY = atof(text.data());
        text = CASCIIUtility::getWord(commandLine, arg++);
        dZ = atof(text.data());
    }
    else
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        lenX = atof(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        lenY = atof(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        lenZ = atof(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        N = atoi(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        dr = atof(text.data());

        useRandomDistribution = true;
    }
    
    // Get list of objects to copy
    do
    {
        text = CASCIIUtility::getWord(commandLine, arg++);
        length = (int)text.size();
        count++;
        
        if(length)
        {
            objectsToCopy.emplace_back(text);
        }
    
    } while(length && (count < 200));

    // Convert objects to proper atom pointers and produce a list: aAtoms to copy,
    if(commandString == "atoms")
    {
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            int iIndex = atoi(objectsToCopy[i].data());
            if((iIndex < 0) || (iIndex >= state_->atoms_.size()))
            {
                printf("Error: found no atom with index %i, aborted copy!", iIndex);
                return;
            }

            CAtom* atom = state_->atoms_[iIndex].get();
            if(!atom)
            {
                printf("Error: found empty atom with index %i, aborted copy!", iIndex);
                return;
            }
        }
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            int iIndex = atoi(objectsToCopy[i].data());
            atomsToCopy.emplace_back(state_->atoms_[iIndex].get());
        }
    }
    else if(commandString == "names")
    {
        std::vector<int> atoms;
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            state_->getAtomsWithID(objectsToCopy[i], atoms);

            for(int j=0; j<atoms.size(); j++)
            {
                atomsToCopy.emplace_back(state_->atoms_[atoms[j]].get());
            }
        }
    }
    else if(commandString == "all")
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            atomsToCopy.emplace_back(state_->atoms_[j].get());
        }
    }
    else if(commandString == "sel")
    {
        for(int j=0; j<state_->atoms_.size(); j++)
        {
            if(state_->atoms_[j]->isSelected())
                atomsToCopy.emplace_back(state_->atoms_[j].get());
        }
    }
    else if(commandString == "resnames")
    {
        std::vector<int> atoms;
        
        for(int i=0; i<objectsToCopy.size(); i++)
        {
            state_->getAtomsWithResname(objectsToCopy[i], atoms);
            
            for(int j=0; j<atoms.size(); j++)
            {
                atomsToCopy.emplace_back(state_->atoms_[atoms[j]].get());
            }
        }
    }
    
    // Copy the atoms, keep only intrinsic atom properties in copy
    if(!useRandomDistribution)
    {
        for(int i=0; i<atomsToCopy.size(); i++)
        {
            CAtom* atomPtr = atomsToCopy[i];

            sX = atomPtr->r_[0].x_;
            sY = atomPtr->r_[0].y_;
            sZ = atomPtr->r_[0].z_;

            CAtom atom;
            atom.copyIntrinsicAtomProperties(*atomPtr);

            for(int iX=0; iX<nX; iX++)
            {
                for(int iY=0; iY<nY; iY++)
                {
                    for(int iZ=0; iZ<nZ; iZ++)
                    {
                        // Make sure we do not make a copy of the atom that is already there
                        if((iX == 0) && (iY == 0) && (iZ == 0)) continue;

                        // Copy to the new position away from pos. of existing atom
                        atom.r_[0].x_ = sX + double(iX)*dX;
                        atom.r_[0].y_ = sY + double(iY)*dY;
                        atom.r_[0].z_ = sZ + double(iZ)*dZ;
                        state_->addAtom(atom);
                    }
                }
            }
        }
    }

    // Perform copy, but use ransom placement of the copies within a given bounding box, and keep only intrinsic atom properties in copy
    else
    {
        // Calculate geometric center
        C3DVector c = CMolTwisterStateTools::getGeometricCenter(atomsToCopy, 0);

        // Find longest vector from c
        double R2_max = 0.0;
        C3DVector R_max;
        for(int i=0; i<atomsToCopy.size(); i++)
        {
            CAtom* atomPtr = atomsToCopy[i];

            C3DVector r = atomPtr->r_[0];
            C3DVector R = r - c;
            double R2 = R.norm2();

            if(R2 > R2_max)
            {
                R2_max = R2;
                R_max = R;
            }
        }

        // Calculate encompassing sphere radius and use this to estimate an appropriate simulation box range to place the centers
        double R_sphere = R_max.norm() + dr;
        double D_sphere = 2.0 * R_sphere;
        double rangeX = lenX - D_sphere;
        double rangeY = lenY - D_sphere;
        double rangeZ = lenZ - D_sphere;

        // Lambda to generate molecular coordinates at r
        std::function<std::shared_ptr<std::vector<C3DVector>>(const C3DVector&, double, double, double, const std::vector<CAtom*>&)> genMolecularCoordinates
                = [c](const C3DVector& r, double alpha, double beta, double gamma, const std::vector<CAtom*>& atomsToCopy)
        {
            C3DVector translate = r - c;
            auto atomsMolecule = std::make_shared<std::vector<C3DVector>>(atomsToCopy.size());
            for(size_t i=0; i<atomsToCopy.size(); i++)
            {
                (*atomsMolecule)[i] = translate + atomsToCopy[i]->r_[0];
                (*atomsMolecule)[i].rotate(r, alpha, beta, gamma);
            }

            return atomsMolecule;
        };

        // Lambda to return non-overlapping molecule if coordinates are not overlapping
        std::function<std::pair<bool, std::shared_ptr<std::vector<C3DVector>>>(const C3DVector&, double, double, double, const std::vector<CAtom*>&, const std::vector<std::shared_ptr<std::vector<C3DVector>>>&)> returnMolIfNonOverlapCoords
                = [dr, genMolecularCoordinates](const C3DVector& r, double alpha, double beta, double gamma, const std::vector<CAtom*>& atomsToCopy, const std::vector<std::shared_ptr<std::vector<C3DVector>>>& placedMolecules)
        {
            const double dr2 = dr*dr;
            auto atomsMolecule1 = genMolecularCoordinates(r, alpha, beta, gamma, atomsToCopy);

            for(const auto& atomsMolecule2 : placedMolecules)
            {
                for(const C3DVector& r1 : *atomsMolecule1)
                {
                    for(const C3DVector& r2 : *atomsMolecule2)
                    {
                        const C3DVector d = r2 - r1;
                        if(d.norm2() < dr2) return std::pair<bool, std::shared_ptr<std::vector<C3DVector>>>(true, nullptr);
                    }
                }
            }

            return std::pair<bool, std::shared_ptr<std::vector<C3DVector>>>(false, atomsMolecule1);
        };

        // Place N copies
        CProgressBar pb;
        std::srand(std::time(nullptr));
        std::vector<std::shared_ptr<std::vector<C3DVector>>> placedMolecules;
        placedMolecules.emplace_back(genMolecularCoordinates(C3DVector(0.0, 0.0, 0.0), 0.0, 0.0, 0.0, atomsToCopy));
        pb.beginProgress("Find non-overlaping centers");
        for(int i=0; i<N; i++)
        {
            // We try to place it max 1000 times
            for(int t=0; t<1000; t++)
            {
                // Try random position and random direction
                const double twoPi = 2.0 * M_PI;
                const double x = double(std::rand()) * rangeX / double(RAND_MAX);
                const double y = double(std::rand()) * rangeY / double(RAND_MAX);
                const double z = double(std::rand()) * rangeZ / double(RAND_MAX);
                const double alpha = double(std::rand()) * twoPi / double(RAND_MAX);
                const double beta = double(std::rand()) * twoPi / double(RAND_MAX);
                const double gamma = double(std::rand()) * twoPi / double(RAND_MAX);

                // Check if we can accept the position and place if OK and move on, else try again
                auto nonOverlappingMolecule = returnMolIfNonOverlapCoords(C3DVector(x, y, z), alpha, beta, gamma, atomsToCopy, placedMolecules);
                if(!nonOverlappingMolecule.first)
                {
                    placedMolecules.emplace_back(nonOverlappingMolecule.second);
                    break;
                }
            }
            pb.updateProgress(i, N);
        }
        pb.endProgress();

        // Delete the first entry of placedMolCenters, so we do not duplicate the original
        placedMolecules.erase(placedMolecules.begin());

        // Create the atomic copies
        pb.beginProgress("Copy atoms");
        for(size_t i=0; i<placedMolecules.size(); i++)
        {
            const auto& molecule = placedMolecules[i];
            if(!molecule) continue;

            for(int j=0; j<atomsToCopy.size(); j++)
            {
                CAtom* atomPtr = atomsToCopy[j];

                CAtom atomCpy;
                atomCpy.copyIntrinsicAtomProperties(*atomPtr);
                atomCpy.r_[0] = (*molecule)[j];
                state_->addAtom(atomCpy);
            }
            pb.updateProgress(i, placedMolecules.size());
        }
        pb.endProgress();

        // Warn if not success and confirm if success
        if(placedMolecules.size() < N)
        {
            printf("Warning: only placed %i of %i molecules, try lowering <dr>!", (int)placedMolecules.size(), N);
        }
        else
        {
            printf("\r\n\tPlaced all %i molecules!\r\n", (int)placedMolecules.size());
        }
    }

    // Update 3D view
    if(state_->view3D_)
    {
        state_->view3D_->requestUpdate(false);
    }
}
