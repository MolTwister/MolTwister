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

#include "CmdMSD.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdMSD::getCmd()
{
    return "msd";
}

std::vector<std::string> CCmdMSD::getCmdLineKeywords()
{
    return { "msd", "resname", "usegeomcent", "numshiftsint0" };
}

std::vector<std::string> CCmdMSD::getCmdHelpLines()
{
    return {
                "msd <DCD filename> <frame from> <frame to> resname <resname> [usegeomcent] [numshiftsint0 <num shifts>] "
           };
}

std::string CCmdMSD::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the mean square displacement (MSD) of the molecules that has the resname <resname>\r\n";
    text+= "\tbetween frames <frame from> to <frame to>. By default, the MSD is calculated based on the\r\n";
    text+= "\tcenter of mass of the tracked molecules (which requires masses of the atoms to be loaded).\r\n";
    text+= "\tHowever, by using the 'usegeomcent' keyword, all masses are set to unity, thus resulting\r\n";
    text+= "\tin the geometric center of the molecules being used as basis for the MSD.\r\n";
    text+= "\r\n";
    text+= "\tTo enable better MSD averaging, the 'numshiftsint0' keyword kan be used to specify <num\r\n";
    text+= "\tshifts> number of shifts to apply to the starting frame, t0. Hence, if <num shifts>=0, no\r\n";
    text+= "\tshifts are applied, if <num shifts>=10, then the MSD is calculated 10 times but with the\r\n";
    text+= "\tstarting frame being moved 0, 1, 2 and up to 10 steps. The average MSD is calculated between\r\n";
    text+= "\tall the shifted calculations.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Num. shifts in t0 = <num shifts>\r\n";
    text+= "\t2. From frame = <frame from>\r\n";
    text+= "\t3. To frame = <frame to>\r\n";
    text+= "\t4. Index MSD\r\n";
    text+= "\t5. <index> <MSD>\r\n";
    text+= "\t6. <index> <MSD>\r\n";
    text+= "\t          .\r\n";
    text+= "\t          .\r\n";
    text+= "\t          .\r\n";
    text+= "\tN+4. <index> <MSD>\r\n";
    text+= "\twhere N is the number of points in the MSD curve.";

    return text;
}

std::string CCmdMSD::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::vector<int>> listOfMolecules;
    CDCDFile dcdFile;
    std::string text;
    bool useGeometricCenter = false;
    int frameFrom, frameTo, numShiftsIn_t0=0;


    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    frameFrom = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    frameTo = atoi(text.data());

    if((frameTo - frameFrom) <= 0)
    {
        lastError_ = "Error: no frames were selected!";
        return lastError_;
    }


    // Find molecules to loop over
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "resname")
    {
        text = CASCIIUtility::getArg(arguments, arg++);

        // Build list of molecules, based on resname, from the atom list. Note! we will not have more molecules
        // than atoms (hence the size of the unorganized list). A molecule is defined based on visible bonds of
        // the 3D view of MolTwister
        std::vector<std::vector<int>> unorganizedListOfMolecules;
        unorganizedListOfMolecules.resize(state_->atoms_.size());

        for(int i=0; i<(int)state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->resname_ == text)
            {
                unorganizedListOfMolecules[state_->atoms_[i]->getMolIndex()].emplace_back(state_->atoms_[i]->getAtomIndex());
            }
        }

        // Pick out the molecules that are available in the unorganized list and build a more organized list
        for(int i=0; i<(int)unorganizedListOfMolecules.size(); i++)
        {
            if(unorganizedListOfMolecules[i].size() == 0) continue;
            listOfMolecules.emplace_back(unorganizedListOfMolecules[i]);
        }
    }
    else
    {
        lastError_ = "Error: expected 'resname'!";
        return lastError_;
    }


    // Get optional switches
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "usegeomcent") useGeometricCenter = true;
    else arg--;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "numshiftsint0")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        numShiftsIn_t0 = atoi(text.data());
        if(numShiftsIn_t0 < 0)
        {
            lastError_ = "Error: number of shifts in t0 cannot be less than zero!";
            return lastError_;
        }
    }
    else arg--;


    // Calculate the MSD (Mean Square Displacement)
    frameTo-= numShiftsIn_t0;
    std::vector<double> msdCurve;
    msdCurve.resize(frameTo - frameFrom, 0.0);
    char progressText[64];
    for(int t0=0; t0<(numShiftsIn_t0+1); t0++)
    {
        sprintf(progressText, "Generating MSD with t0=%i (%i/%i)", frameFrom + t0 + 1, t0 + 1, numShiftsIn_t0 + 1);
        pb.beginProgress(progressText);

        std::vector<C3DVector> aRcAtTime0;
        aRcAtTime0.resize(listOfMolecules.size());
        for(int t=frameFrom+t0; t<frameTo+t0; t++)
        {
            dcdFile.gotoRecord(t);
            int maxNumCoordinates = dcdFile.getNumCoordinatesInRecord();

            double sum = 0.0;
            for(int i=0; i<(int)listOfMolecules.size(); i++)
            {
                // Get center of mass ($\mathbf{R}_c = \frac{1}{M_{tot}} \sum_{n=1}^{N_{atoms}} m_n \mathbf{r}_n$)
                // or geomeric center (where $m_n\equiv 1$ and $M = N_{atoms}$). $m_n$ is the mass of atom $n$.
                double M = 0.0;
                C3DVector Rc;
                for(int n=0; n<(int)listOfMolecules[i].size(); n++)
                {
                    double m = 1.0;
                    if(!useGeometricCenter) m = state_->atoms_[listOfMolecules[i][n]]->m_;
                    if(fabs(m) < 1E-5)
                    {
                        lastError_ = "Error: found zero mass (use 'usegeomcent' or load masses)!";
                        return lastError_;
                    }
                    M+= m;

                    if(listOfMolecules[i][n] >= maxNumCoordinates)
                    {
                        printf("Error: coordinate index %i does not exist in DCD file (num coordinates in DCD = %i)!", listOfMolecules[i][n], maxNumCoordinates);
                        lastError_ = std::string("Error: coordinate index ") + std::to_string(listOfMolecules[i][n]) + std::string(" does not exist in DCD file (num coordinates in DCD = ") + std::to_string(maxNumCoordinates) + std::string(")!");
                        return lastError_;
                    }
                    Rc+= dcdFile.getCoordinate(listOfMolecules[i][n])*m;
                }
                if(M != 0.0) Rc*= (1.0 / M);

                // Store $R_c(t=0)$ for all molecules
                if(t == (frameFrom+t0)) aRcAtTime0[i] = Rc;

                // Calculate MSD contribution for molecule i
                C3DVector deltaRc = Rc - aRcAtTime0[i];
                sum+= deltaRc.norm2();
            }

            if(listOfMolecules.size() == 0)
            {
                lastError_ = std::string("Error: no molecules with specified resname was found in frame ") + std::to_string(t) + std::string("!");
                return lastError_;
            }

            // Add to the MSD curve in the correct bin
            sum/= double(listOfMolecules.size());
            msdCurve[t-frameFrom-t0]+= sum;

            pb.updateProgress(t-frameFrom-t0, frameTo - frameFrom);
        }

        pb.endProgress();
    }

    // Average each bin in the MSD curve over the number of shifts done in t0
    for(int i=0; i<(int)msdCurve.size(); i++)
    {
        msdCurve[i]/= double(numShiftsIn_t0 + 1);
    }

    // Plot the results
    fprintf(stdOut_, "\tNum. shifts in t0 = %i\r\n", numShiftsIn_t0);
    fprintf(stdOut_, "\tFrom frame = %i\r\n", frameFrom);
    fprintf(stdOut_, "\tTo frame = %i\r\n\r\n", frameTo+numShiftsIn_t0);
    fprintf(stdOut_, "\t%-15s%-15s\r\n", "Index", "MSD");
    for(int i=0; i<(int)msdCurve.size(); i++)
    {
        fprintf(stdOut_, "\t%-15i% -15.8f\r\n", i, msdCurve[i]);
    }

    return lastError_;
}
