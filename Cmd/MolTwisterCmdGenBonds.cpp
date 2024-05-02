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
#include "Utilities/3DRect.h"
#include "MolTwisterCmdGenBonds.h"
#include "Tools/MolTwisterStateTools.h"

void CCmdGenBonds::onAddKeywords()
{
    addKeyword("genbonds");
    addKeyword("verbose");
    addKeyword("nonverbose");
    addKeyword("pbcdetect");
    addKeyword("atomicunwrap");
    addKeyword("minr");
    addKeyword("ignore");
}

std::string CCmdGenBonds::getHelpString() const
{ 
    std::string  text;
    
    text+= "\tUsage: genbonds [minr <minimum R>] [pbcdetect] [verbose] [atomicunwrap]\r\n";
    text+= "\r\n";
    text+= "\tGenerate bonds between atoms. A bond is defined as any bond length R that\r\n";
    text+= "\tsatisfies 'minimum R' < r < r1+r2+0.4. In addition, the number of bonds that\r\n";
    text+= "\tare connected to C, N, P, S atoms are restricted to 4. The default minimum\r\n";
    text+= "\tbond length limit is 0.8AA. The parmeters inside square brackets are optional.\r\n";
    text+= "\tIf the verbose keyword is used, then each detected bond will be displayed as\r\n";
    text+= "\tthey are detected.\r\n";
    text+= "\r\n";
    text+= "\tUsing the pbcdetect keyword results in bonds across periodic images defined\r\n";
    text+= "\tby the current PBCs (i.e. Periodic Boundary Conditions). 'atomicunwrap'\r\n";
    text+= "\tresults in molecules diveded by the PBCs being wrapped to one side of the\r\n";
    text+= "\tperiodic boundary that divide the molecule.\r\n";
    text+= "\r\n";
    text+= "\tTo ignore bonds to atoms use 'ignore' followed by a comma separated list of\r\n";
    text+= "\tatoms (without any space).";
    
    return text;
}

void CCmdGenBonds::execute(std::string commandLine)
{
    std::vector<std::string> bondAtomsToIgnore;
    C3DRect* pbcPtr = nullptr;
    C3DRect pbc;
    std::string text;
    double bondLength = 0.8;
    bool verbose = false;
    bool applyAtomicUnwrap = false;
    
    if(!state_) return;
    if(!state_->view3D_) return;
    
    pbc = state_->view3D_->calcPBC();

    for(int i=1; i<=20; i++)
    {
        text = CASCIIUtility::getWord(commandLine, i);
        CASCIIUtility::removeWhiteSpace(text);
        if(text == "verbose") verbose = true;
        if(text == "pbcdetect") pbcPtr = &pbc;
        if(text == "atomicunwrap") { pbcPtr = &pbc; applyAtomicUnwrap = true; }
        if(text == "minr")
        {
            text = CASCIIUtility::getWord(commandLine, ++i);
            CASCIIUtility::removeWhiteSpace(text);
            bondLength = atof(text.data());
        }
        if(text == "ignore")
        {
            text = CASCIIUtility::getWord(commandLine, ++i);
            CASCIIUtility::removeWhiteSpace(text);
            bondAtomsToIgnore = CASCIIUtility::getWords(text, ",");
        }
    }
    
    fprintf(stdOut_, "\r\n\tBond lenght criteria 1: R in [%.4f, r1+r2+0.4]\r\n", bondLength);
    
    CMolTwisterStateTools(state_, stdOut_).generateBonds(bondLength, verbose, true, -1, pbcPtr, &bondAtomsToIgnore);
    if(applyAtomicUnwrap) CMolTwisterStateTools(state_, stdOut_).atomicUnwrap(pbc);
    state_->view3D_->requestUpdate(false);
}
