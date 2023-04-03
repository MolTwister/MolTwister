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

#include "CmdQBal.h"
#include "../../Utilities/ASCIIUtility.h"
#include <math.h>

std::string CCmdQBal::getCmd()
{
    return "qbal";
}

std::vector<std::string> CCmdQBal::getCmdLineKeywords()
{
    return { "qbal" };
}

std::vector<std::string> CCmdQBal::getCmdHelpLines()
{
    return {
                "qbal <group of atoms to modify>"
           };
}

std::string CCmdQBal::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tGiven a loaded dataset, this command will calculate scaling factors for the atoms within the given set of\r\n";
    text+= "\tatoms to modify. The <group of atoms to modify> parameter can be one of the following:\r\n";
    text+= "\t* sel - the group of atoms to modify is the current visual selection\r\n";
    text+= "\tThe calculated scaling factors for the charges are such that, if applied to the corresponding charges, the\r\n";
    text+= "\ttotal charge of the loaded system is charge balanced. Moreover, the scaling factors are minimized, such that\r\n";
    text+= "\tthey are as small as possible (given the usual L2-norm).\r\n";
    text+= "\r\n";
    text+= "\tMore detailed information about calculating the scaling factors is given in Refs. 1 and 2 (see below).\r\n";
    text+= "\t\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Charge 'Multiply with' 'Old charge' 'New charge' 'Change in %'\r\n";
    text+= "\t2. <charge type> <correction factor> <old charge value> <charge value after multiplication with correction factor> <change in % from old charge>\r\n";
    text+= "\t3. <charge type> <correction factor> <old charge value> <charge value after multiplication with correction factor> <change in % from old charge>\r\n";
    text+= "\t         .\r\n";
    text+= "\t         .\r\n";
    text+= "\t         .\r\n";
    text+= "\tN+1. <charge type> <correction factor> <old charge value> <charge value after multiplication with correction factor> <change in % from old charge>\r\n";
    text+= "\twhere N is the number of charge types (e.g., H, O, C8) within the group of atoms to modify.\r\n";
    text+= "\r\n";
    text+= "\r\n";
    text+= "\t[1] Olsen, R. and Kvamme, B. (2019) ‘Effects of glycol on adsorption dynamics of idealized water droplets on LTA‐3A zeolite surfaces’, AIChE Journal, 65(5), p. e16567. doi: 10.1002/aic.16567.\r\n";
    text+= "\t[2] Olsen, R. et al. (2016) ‘Effects of Sodium Chloride on Acidic Nanoscale Pores Between Steel and Cement’, The Journal of Physical Chemistry C, 120(51), pp. 29264–29271. doi: 10.1021/acs.jpcc.6b10043.";

    return text;
}

std::string CCmdQBal::execute(std::vector<std::string> arguments)
{
    lastError_ = "";
    size_t arg = 0;

    // Get the type of atom filtering
    std::string text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "sel")
    {
        std::vector<CAtom*> selAtoms;
        std::vector<std::string> types;
        std::vector<double> q_ixn_i, q_ixn_i2;
        std::vector<double> firstQ;
        double Q = 0.0;
        double lambda = 0.0;
        bool s1 = false;
        bool s2 = false;

        // Find total charge, Q
        for(int i=0; i<state_->atoms_.size(); i++)
            Q+= state_->atoms_[i]->Q_;

        // Get types of selected atoms, q_i*n_i for each type, and first charge of each type
        findAtomTypesInSel(types);
        for(int i=0; i<types.size(); i++)
        {
            // Find instances of this type
            findInstancesOfAtomType(types[i], true, selAtoms);

            // Store the first charge from all instances
            if(selAtoms.size() > 0) firstQ.emplace_back(selAtoms[0]->Q_);
            else                    firstQ.emplace_back(0.0);

            // Store the total charge of the selection of this type
            double q = 0.0;
            for(int j=0; j<selAtoms.size(); j++)
                q+= selAtoms[j]->Q_;

            q_ixn_i.emplace_back(q);
            q_ixn_i2.emplace_back(q*q);
        }

        // Newton-Raphson to find lambda in case of non-negative scaling
        for(int I=0; I<100; I++)
        {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for(int i=0; i<q_ixn_i.size(); i++)
            {
                double g = 1.0 - lambda*q_ixn_i[i]/2.0;
                sum1+= (q_ixn_i[i] * (1.0 - fabs(g)));
                sum2+= q_ixn_i2[i] * sgn(g);
            }

            double f = Q - sum1;
            if(fabs(f) < 1.0E-12)
            {
                s1 = true;
                break;
            }
            double Df = -sum2 / 2.0;
            if(Df == 0.0) break;
            lambda = lambda - (f / Df);
        }

        // If we failed to find non-negative scaling then...
        if(!s1)
        {
            double D = 0.0;
            for(int i=0; i<q_ixn_i.size(); i++)
            {
                D+= q_ixn_i2[i];
            }

            if(D == 0.0)
            {
                lastError_ = "Error: could not balance the charges, division by zero!";
                return lastError_;
            }

            lambda = 2.0 * Q / D;
            s2 = true;
        }

        fprintf(stdOut_, "\r\n\t%-10s%-20s%-20s%-15s%-15s\r\n\t----------------------------------------------------------------------------\r\n", "Charge", "Multiply with", "Old charge", "New charge", "Change in %");

        for(int i=0; i<types.size(); i++)
        {
            double alpha_i = 1.0 - lambda * q_ixn_i[i] / 2.0;

            if(s1) alpha_i = fabs(alpha_i);
            fprintf(stdOut_, "\t%-10s% -20.8f% -20.8f% -15.8f% -15.2f\r\n", types[i].data(), alpha_i, firstQ[i], alpha_i*firstQ[i], fabs((alpha_i - 1.0))*100.0);
        }

        fprintf(stdOut_, "\r\n\tLp-norm = %i\r\n", (int)2.0);
        if(s2) fprintf(stdOut_, "\r\n\tWarning: could not find all non-negative charge scale factors!\r\n");
    }

    return lastError_;
}

void CCmdQBal::findAtomTypesInSel(std::vector<std::string>& types) const
{
    types.clear();

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(!state_->atoms_[i]->isSelected()) continue;

        // Get ID
        std::string ID =  state_->atoms_[i]->getID();

        // Check if we already added it
        bool inTypes = false;
        for(int j=0; j<types.size(); j++)
        {
            if(ID == types[j]) inTypes = true;
        }

        // If not add it
        if(!inTypes)
        {
            types.emplace_back(ID);
        }
    }
}

void CCmdQBal::findInstancesOfAtomType(std::string type, bool includeSelOnly, std::vector<CAtom*>& atoms) const
{
    atoms.clear();

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(includeSelOnly && !state_->atoms_[i]->isSelected()) continue;

        std::string ID = state_->atoms_[i]->getID();
        if(ID == type)
        {
            atoms.emplace_back(state_->atoms_[i].get());
        }
    }
}
