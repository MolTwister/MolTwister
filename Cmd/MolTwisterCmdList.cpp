//
// Copyright (C) 2021 Richard Olsen.
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
#include "MolTwisterCmdList.h"

void CCmdList::onAddKeywords()
{
    addKeyword("list");
    addKeyword("all");
    addKeyword("mol");
    addKeyword("ff");
    addKeyword("latexff");
    addKeyword("longtable");
    addKeyword("longtable");
}

std::string CCmdList::getHelpString() const
{ 
    std::string text;
    
    text+= "\tUsage: list <filter>\r\n";
    text+= "\r\n";
    text+= "\tLists atomic properties based on the specified filter. The allowed filters\r\n";
    text+= "\tare:\r\n";
    text+= "\r\n";
    text+= "\t       * all                    :   List all\r\n";
    text+= "\t       * mol <N>                :   List atoms in molecule N\r\n";
    text+= "\t       * ff                     :   List only force-field parameters\r\n";
    text+= "\t       * latexff [longtable]    :   Make LaTeX force-field tables. 'longtable'\r\n";
    text+= "\t                                    creates tables across several pages\r\n";
    
    return text;
}

void CCmdList::execute(std::string commandLine)
{
    int arg = 1;
    std::string text;
     
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "all")
    {
        parseAllCommand(commandLine, arg);
    }
    else if(text == "mol")
    {
        parseMolCommand(commandLine, arg);
    }
    else if(text == "ff")
    {
        parseFFCommand(commandLine, arg);
    }
    else if(text == "latexff")
    {
        parseLatexCommand(commandLine, arg);
    }
    else
    {
        printf("Syntax Error: Second argument should be the type of list!");
    }
}

void CCmdList::parseAllCommand(std::string, int&)
{
    printListHeaderAtoms();
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        printListAtomicItem(state_->atoms_[i].get());
    }

    if(state_->mdFFNonBondedList_.size() > 0)
    {
        fprintf(stdOut_, "\r\n");
        printListHeaderMDNonBond();
        for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
        {
            printListNonBondedItem(state_->mdFFNonBondedList_.get(i), i);
        }
    }

    if(state_->mdFFBondList_.size() > 0)
    {
        fprintf(stdOut_, "\r\n");
        printListHeaderMDBonds();
        for(int i=0; i<state_->mdFFBondList_.size(); i++)
        {
            printListBondItem(state_->mdFFBondList_.get(i), i);
        }
    }

    if(state_->mdFFAngleList_.size() > 0)
    {
        fprintf(stdOut_, "\r\n");
        printListHeaderMDAngles();
        for(int i=0; i<state_->mdFFAngleList_.size(); i++)
        {
            printListAngleItem(state_->mdFFAngleList_.get(i), i);
        }
    }

    if(state_->mdFFDihList_.size() > 0)
    {
        fprintf(stdOut_, "\r\n");
        printListHeaderMDDihedrals();
        for(int i=0; i<state_->mdFFDihList_.size(); i++)
        {
            printListDihItem(state_->mdFFDihList_.get(i), i);
        }
    }
}

void CCmdList::parseMolCommand(std::string commandLine, int& arg)
{
    int molIndex;

    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    molIndex = atoi(text.data());

    printListHeaderAtoms();
    for(int i=0; i<(int)state_->atoms_.size(); i++)
    {
        if(molIndex == state_->atoms_[i]->getMolIndex())
            printListAtomicItem(state_->atoms_[i].get());
    }
}

void CCmdList::parseFFCommand(std::string, int&)
{
    bool firstPrinted = false;

    if(state_->mdFFNonBondedList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printListHeaderMDNonBond();
        for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
        {
            printListNonBondedItem(state_->mdFFNonBondedList_.get(i), i);
        }

        firstPrinted = true;
    }

    if(state_->mdFFBondList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printListHeaderMDBonds();
        for(int i=0; i<state_->mdFFBondList_.size(); i++)
        {
            printListBondItem(state_->mdFFBondList_.get(i), i);
        }

        firstPrinted = true;
    }

    if(state_->mdFFAngleList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printListHeaderMDAngles();
        for(int i=0; i<state_->mdFFAngleList_.size(); i++)
        {
            printListAngleItem(state_->mdFFAngleList_.get(i), i);
        }

        firstPrinted = true;
    }

    if(state_->mdFFDihList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printListHeaderMDDihedrals();
        for(int i=0; i<state_->mdFFDihList_.size(); i++)
        {
            printListDihItem(state_->mdFFDihList_.get(i), i);
        }

        firstPrinted = true;
    }
}

void CCmdList::parseLatexCommand(std::string commandLine, int& arg)
{
    bool firstPrinted = false;
    bool useLongTable = false;

    std::string text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "longtable") useLongTable = true;

    if(state_->mdFFNonBondedList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printLaTeXListHeaderMDNonBond(useLongTable);
        for(int i=0; i<state_->mdFFNonBondedList_.size(); i++)
        {
            printLaTeXListNonBondedItem(state_->mdFFNonBondedList_.get(i), i, useLongTable);
        }
        printLaTeXListFooterMDNonBond(useLongTable);

        firstPrinted = true;
    }

    if(state_->mdFFBondList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printLaTeXListHeaderMDBonds(useLongTable);
        for(int i=0; i<state_->mdFFBondList_.size(); i++)
        {
            printLaTeXListBondItem(state_->mdFFBondList_.get(i), i, useLongTable);
        }
        printLaTeXListFooterMDBonds(useLongTable);

        firstPrinted = true;
    }

    if(state_->mdFFAngleList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printLaTeXListHeaderMDAngles(useLongTable);
        for(int i=0; i<state_->mdFFAngleList_.size(); i++)
        {
            printLaTeXListAngleItem(state_->mdFFAngleList_.get(i), i, useLongTable);
        }
        printLaTeXListFooterMDAngles(useLongTable);

        firstPrinted = true;
    }

    if(state_->mdFFDihList_.size() > 0)
    {
        if(firstPrinted) fprintf(stdOut_, "\r\n");
        printLaTeXListHeaderMDDihedrals(useLongTable);
        for(int i=0; i<state_->mdFFDihList_.size(); i++)
        {
            printLaTeXListDihItem(state_->mdFFDihList_.get(i), i, useLongTable);
        }
        printLaTeXListFooterMDDihedrals(useLongTable);

        firstPrinted = true;
    }

    if(firstPrinted) fprintf(stdOut_, "\r\n");
    printLaTeXListCharges(useLongTable);
}

void CCmdList::printLaTeXListCharges(bool useLongTable) const
{
    std::vector<std::string> listOfAtomTypes;
    std::vector<std::string> listOfResnames;
    CAtom* atomPtr;
    
    state_->searchForAtomTypes(listOfAtomTypes, &listOfResnames);
    
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\\begin{longtable}{l r l r l r}\r\n");
        fprintf(stdOut_, "\t\t\\toprule Atom & $Q$ & Atom & $Q$ & Atom & $Q$ \\\\ \\midrule\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\\begin{table}[t!]\r\n");
        fprintf(stdOut_, "\t\t\\centering\r\n");
        fprintf(stdOut_, "\t\t\t\\begin{tabular}{l r l r l r}\r\n");
        fprintf(stdOut_, "\t\t\t\t\\toprule Atom & $Q$ & Atom & $Q$ & Atom & $Q$ \\\\ \\midrule\r\n");
    }
    
    for(int i=0; i<(int)listOfAtomTypes.size(); i++)
    {
        if(useLongTable)   fprintf(stdOut_, "\t\t\t");
        else               fprintf(stdOut_, "\t\t\t\t\t");
        
        for(int j=0; j<3; j++)
        {
            if(j != 0) fprintf(stdOut_, "& ");
            if(i < (int)listOfAtomTypes.size())
            {
                atomPtr = state_->getFirstOccurenceOf(listOfAtomTypes[i]);
                fprintf(stdOut_, "%s & %.8f ", listOfAtomTypes[i].data(), atomPtr ? atomPtr->Q_ : 0.0);
            }
            else
            {
                fprintf(stdOut_, " &  ");
            }
            
            i++;
        }
        i--;
        
        fprintf(stdOut_, "\\\\ \r\n");
    }

    if(useLongTable)
    {
        fprintf(stdOut_, "\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\\caption{Partial charges assigned to each atom type.}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{longtable}\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\t\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\t\\end{tabular}\r\n");
        fprintf(stdOut_, "\t\t\\caption{Partial charges assigned to each atom type.\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{table}\r\n");
    }
}

void CCmdList::printListHeaderAtoms() const
{
    fprintf(stdOut_, "\r\n\r\n\t-------------------------------------------------------------- ATOMS ------------------------------------------------------------------------------------------------\r\n");
    fprintf(stdOut_, "\t%-10s%-11s%-9s%-9s%-15s%-15s%-15s%-15s%-15s%-7s%-15s%-10s%-28s\r\n", "At.Index", "Mol.Index", "At.Name", "resname", " x", " y", " z", "Q", "m", "mobile", "sigma (vdW)", "Selected", "Bonded to:[At.Index,At.ID]");
    fprintf(stdOut_, "\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\r\n");
}

void CCmdList::printListHeaderMDNonBond() const
{
    fprintf(stdOut_, "\r\n\r\n\t------------------------------------------------------------ NON-BONDED ---------------------------------------------------------------------------------------------\r\n");
    fprintf(stdOut_, "\t%-10s%-10s%-10s%-15s%-20s%-30s\r\n", "Index", "Atom 1", "Atom 2", "Type", "Comments", "Arguments");
    fprintf(stdOut_, "\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\r\n");
}

void CCmdList::printLaTeXListHeaderMDNonBond(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\\begin{longtable}{l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\\toprule Atom $i$ & Atom $j$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\\begin{table}[t!]\r\n");
        fprintf(stdOut_, "\t\t\\centering\r\n");
        fprintf(stdOut_, "\t\t\t\\begin{tabular}{l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\t\t\\toprule Atom $i$ & Atom $j$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
}

void CCmdList::printLaTeXListFooterMDNonBond(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{longtable}\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\t\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\t\\end{tabular}\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{table}\r\n");
    }
}

void CCmdList::printListHeaderMDBonds() const
{
    fprintf(stdOut_, "\r\n\r\n\t-------------------------------------------------------------- BONDS ------------------------------------------------------------------------------------------------\r\n");
    fprintf(stdOut_, "\t%-10s%-10s%-10s%-15s%-20s%-20s%-30s\r\n", "Index", "Atom 1", "Atom 2", "Type", "Bond crit.", "Comments", "Arguments");
    fprintf(stdOut_, "\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\r\n");
}

void CCmdList::printLaTeXListHeaderMDBonds(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\\begin{longtable}{l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\\toprule Atom $i$ & Atom $j$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\\begin{table}[t!]\r\n");
        fprintf(stdOut_, "\t\t\\centering\r\n");
        fprintf(stdOut_, "\t\t\t\\begin{tabular}{l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\t\t\\toprule Atom $i$ & Atom $j$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
}

void CCmdList::printLaTeXListFooterMDBonds(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{longtable}\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\t\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\t\\end{tabular}\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{table}\r\n");
    }
}

void CCmdList::printListHeaderMDAngles() const
{
    fprintf(stdOut_, "\r\n\r\n\t-------------------------------------------------------------- ANGLES -----------------------------------------------------------------------------------------------\r\n");
    fprintf(stdOut_, "\t%-10s%-10s%-10s%-10s%-15s%-20s%-20s%-30s\r\n", "Index", "Atom 1", "Atom 2", "Atom 3", "Type", "Bond crit.", "Comments", "Arguments");
    fprintf(stdOut_, "\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\r\n");
}

void CCmdList::printLaTeXListHeaderMDAngles(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\\begin{longtable}{l l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\\toprule Atom $i$ & Atom $j$ & Atom $k$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\\begin{table}[t!]\r\n");
        fprintf(stdOut_, "\t\t\\centering\r\n");
        fprintf(stdOut_, "\t\t\t\\begin{tabular}{l l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\t\t\\toprule Atom $i$ & Atom $j$ & Atom $k$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
}

void CCmdList::printLaTeXListFooterMDAngles(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{longtable}\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\t\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\t\\end{tabular}\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{table}\r\n");
    }
}

void CCmdList::printListHeaderMDDihedrals() const
{
    fprintf(stdOut_, "\r\n\r\n\t------------------------------------------------------------- DIHEDRALS ---------------------------------------------------------------------------------------------\r\n");
    fprintf(stdOut_, "\t%-10s%-10s%-10s%-10s%-10s%-15s%-20s%-20s%-30s\r\n", "Index", "Atom 1", "Atom 2", "Atom 3", "Atom 4", "Type", "Bond crit.", "Comments", "Arguments");
    fprintf(stdOut_, "\t---------------------------------------------------------------------------------------------------------------------------------------------------------------------\r\n");
}

void CCmdList::printLaTeXListHeaderMDDihedrals(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\\begin{longtable}{l l l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\\toprule Atom $i$ & Atom $j$ & Atom $k$ & Atom $l$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\\begin{table}[t!]\r\n");
        fprintf(stdOut_, "\t\t\\centering\r\n");
        fprintf(stdOut_, "\t\t\t\\begin{tabular}{l l l l l l l}\r\n");
        fprintf(stdOut_, "\t\t\t\t\\toprule Atom $i$ & Atom $j$ & Atom $k$ & Atom $l$ & Type & Parameters & Comments \\\\ \\midrule\r\n");
    }
}

void CCmdList::printLaTeXListFooterMDDihedrals(bool useLongTable) const
{
    if(useLongTable)
    {
        fprintf(stdOut_, "\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{longtable}\r\n");
    }
    else
    {
        fprintf(stdOut_, "\t\t\t\t\\bottomrule\r\n");
        fprintf(stdOut_, "\t\t\t\\end{tabular}\r\n");
        fprintf(stdOut_, "\t\t\\caption{\\textbf{TODO:} place caption text here. REMEMBER TO SPECIFY UNITS!}\r\n");
        fprintf(stdOut_, "\t\t\\label{tabN}\r\n");
        fprintf(stdOut_, "\t\\end{table}\r\n");
    }
}

void CCmdList::printListAtomicItem(const CAtom* atom) const
{
    C3DVector r;
       
    std::string ID1 = atom->getID();
    if(state_->currentFrame_ < (int)atom->r_.size()) r = atom->r_[state_->currentFrame_];
    else if(atom->r_.size() > 0)                r = atom->r_[0];
    else                                        return;
    
    fprintf(stdOut_, "\t%-10i%-11i%-9s%-9s% -15.4f% -15.4f% -15.4f% -15.4f% -15.4f%-7s% -15.4f%-10s\t", atom->getAtomIndex(), atom->getMolIndex(), ID1.data(),
            atom->resname_.data(), r.x_, r.y_, r.z_, atom->Q_, atom->m_, atom->isMobile_ ? " " : "rigid", atom->sigma_, atom->isSelected() ? "[x]" : "[ ]");
    
    for(int j=0; j<atom->getNumBonds(); j++)
    {
        std::string ID2 = atom->getBondDest(j)->getID();
        fprintf(stdOut_, "[%i,%s]", atom->getBondDest(j)->getAtomIndex(), ID2.data());
    }
    
    fprintf(stdOut_, "\r\n");
}

void CCmdList::printListNonBondedItem(const CMDFFNonBonded* pNonBonded, int index) const
{
    fprintf(stdOut_, "\t%-10i%-10s%-10s%-15s%-20s%-30s\t", index, pNonBonded->getAtomInBond(0).data(), pNonBonded->getAtomInBond(1).data(),
            pNonBonded->getFFType().data(), pNonBonded->getComments().data(), pNonBonded->getArguments(false).data());
    fprintf(stdOut_, "\r\n");
}

void CCmdList::printLaTeXListNonBondedItem(const CMDFFNonBonded* nonBonded, int, bool useLongTable) const
{
    if(!useLongTable) fprintf(stdOut_, "\t\t");
    fprintf(stdOut_, "\t\t\t%s & %s & %s & %s & %s \\\\ \r\n", nonBonded->getAtomInBond(0).data(), nonBonded->getAtomInBond(1).data(),
            nonBonded->getFFType().data(), nonBonded->getArgumentsLaTeX(false).data(), nonBonded->getComments().data());
}

void CCmdList::printListBondItem(const CMDFFBond* bond, int index) const
{
    char bondScan[256];
    double R0;
    
    if(bond->getBondDetectionCriteria(R0) == CMDFFBond::critAllRLessThan_R0) sprintf(bondScan, "all r<%.4f", R0);
    else if(bond->getBondDetectionCriteria(R0) == CMDFFBond::critMolRLessThan_R0) sprintf(bondScan, "mol r<%.4f", R0);
    else if(bond->getBondDetectionCriteria(R0) == CMDFFBond::critOnlyVisibleBonds) sprintf(bondScan, "only visible");
    else if(bond->getBondDetectionCriteria(R0) == CMDFFBond::critOnly14Bonds) sprintf(bondScan, "only 1-4");
    else sprintf(bondScan, "unknown");

    fprintf(stdOut_, "\t%-10i%-10s%-10s%-15s%-20s%-20s%-30s\t", index, bond->getAtomInBond(0).data(), bond->getAtomInBond(1).data(),
            bond->getFFType().data(), bondScan, bond->getComments().data(), bond->getArguments(false).data());
    fprintf(stdOut_, "\r\n");
}

void CCmdList::printLaTeXListBondItem(const CMDFFBond* bond, int, bool useLongTable) const
{
    if(!useLongTable) fprintf(stdOut_, "\t\t");
    fprintf(stdOut_, "\t\t\t%s & %s & %s & %s & %s \\\\ \r\n", bond->getAtomInBond(0).data(), bond->getAtomInBond(1).data(),
            bond->getFFType().data(), bond->getArgumentsLaTeX(false).data(), bond->getComments().data());
}

void CCmdList::printListAngleItem(const CMDFFAngle* angle, int index) const
{
    char bondScan[256];
    double R0;
    
    if(angle->getBondDetectionCriteria(R0) == CMDFFAngle::critAllRLessThan_R0) sprintf(bondScan, "all r<%.4f", R0);
    else if(angle->getBondDetectionCriteria(R0) == CMDFFAngle::critMolRLessThan_R0) sprintf(bondScan, "mol r<%.4f", R0);
    else if(angle->getBondDetectionCriteria(R0) == CMDFFAngle::critOnlyVisibleBonds) sprintf(bondScan, "only visible");
    else sprintf(bondScan, "unknown");
    
    fprintf(stdOut_, "\t%-10i%-10s%-10s%-10s%-15s%-20s%-20s%-30s\t", index, angle->getAtomInBond(0).data(), angle->getAtomInBond(1).data(),
            angle->getAtomInBond(2).data(), angle->getFFType().data(), bondScan, angle->getComments().data(), angle->getArguments(false).data());
    fprintf(stdOut_, "\r\n");
}

void CCmdList::printLaTeXListAngleItem(const CMDFFAngle* angle, int, bool useLongTable) const
{
    if(!useLongTable) fprintf(stdOut_, "\t\t");
    fprintf(stdOut_, "\t\t\t%s & %s & %s & %s & %s & %s \\\\ \r\n", angle->getAtomInBond(0).data(), angle->getAtomInBond(1).data(), angle->getAtomInBond(2).data(),
            angle->getFFType().data(), angle->getArgumentsLaTeX(false).data(), angle->getComments().data());
}

void CCmdList::printListDihItem(const CMDFFDih* dih, int index) const
{
    char bondScan[256];
    double R0;
    
    if(dih->getBondDetectionCriteria(R0) == CMDFFDih::critAllRLessThan_R0) sprintf(bondScan, "all r<%.4f", R0);
    else if(dih->getBondDetectionCriteria(R0) == CMDFFDih::critMolRLessThan_R0) sprintf(bondScan, "mol r<%.4f", R0);
    else if(dih->getBondDetectionCriteria(R0) == CMDFFDih::critOnlyVisibleBonds) sprintf(bondScan, "only visible");
    else sprintf(bondScan, "unknown");

    fprintf(stdOut_, "\t%-10i%-10s%-10s%-10s%-10s%-15s%-20s%-20s%-30s\t", index, dih->getAtomInBond(0).data(), dih->getAtomInBond(1).data(), dih->getAtomInBond(2).data(), dih->getAtomInBond(3).data(),
            dih->getFFType().data(), bondScan, dih->getComments().data(), dih->getArguments(false).data());
    fprintf(stdOut_, "\r\n");
}

void CCmdList::printLaTeXListDihItem(const CMDFFDih* dih, int, bool useLongTable) const
{
    if(!useLongTable) fprintf(stdOut_, "\t\t");
    fprintf(stdOut_, "\t\t\t%s & %s & %s & %s & %s & %s & %s \\\\ \r\n", dih->getAtomInBond(0).data(), dih->getAtomInBond(1).data(), dih->getAtomInBond(2).data(), dih->getAtomInBond(3).data(),
            dih->getFFType().data(), dih->getArgumentsLaTeX(false).data(), dih->getComments().data());
}
