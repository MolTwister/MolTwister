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

#include "ParsingTools.h"
#include "../../Utilities/ASCIIUtility.h"

bool CParsingTools::retrieve1to4BondSepCoeffs(std::vector<std::string> arguments, size_t& arg, double* a1to4BondSepCoeffs) const
{
    std::string text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "1to4coeffs")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        a1to4BondSepCoeffs[0] = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        a1to4BondSepCoeffs[1] = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        a1to4BondSepCoeffs[2] = atof(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        a1to4BondSepCoeffs[3] = atof(text.data());
    }
    else
    {
        arg--;
        return false;
    }

    return true;
}

std::pair<bool, std::string> CParsingTools::retrieveDihedralAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, CAtom** atom4Ptr, int* atomIndices) const
{
    std::string text;
    std::string lastError;
    bool foundIndex = false;
    bool succeeded = false;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[0] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[1] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[2] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[3] = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;

        text = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeDihedral))
        {
            CVarDihedral* p = (CVarDihedral*)varPtr;

            atomIndices[0] = p->atomIndex1_;
            atomIndices[1] = p->atomIndex2_;
            atomIndices[2] = p->atomIndex3_;
            atomIndices[3] = p->atomIndex4_;
            foundIndex = true;
        }
        else
        {
            lastError = "Variable missing or not of type 'dihedral'!";
        }
    }
    else
    {
        lastError = "Syntax Error: Third argument should be 'id' or 'var'!";
    }


    if(foundIndex)
    {
        if((atomIndices[0] < state_->atoms_.size()) && (atomIndices[1] < state_->atoms_.size())
           && (atomIndices[2] < state_->atoms_.size()) && (atomIndices[3] < state_->atoms_.size()))
        {
            *atom1Ptr = state_->atoms_[atomIndices[0]].get();
            *atom2Ptr = state_->atoms_[atomIndices[1]].get();
            *atom3Ptr = state_->atoms_[atomIndices[2]].get();
            *atom4Ptr = state_->atoms_[atomIndices[3]].get();

            succeeded = true;
        }
        else
        {
            lastError = "Invalid atom indices!";
        }
    }

    return std::pair<bool, std::string>(succeeded, lastError);
}

std::pair<bool, std::string> CParsingTools::retrieveAngleAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, CAtom** atom3Ptr, int* atomIndices) const
{
    std::string lastError;
    std::string text;
    bool foundIndex = false;
    bool succeeded = false;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[0] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[1] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[2] = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;

        text = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeAngle))
        {
            CVarAngle* p = (CVarAngle*)varPtr;

            atomIndices[0] = p->atomIndex1_;
            atomIndices[1] = p->atomIndex2_;
            atomIndices[2] = p->atomIndex3_;
            foundIndex = true;
        }
        else
        {
            lastError = "Variable missing or not of type 'angle'!";
        }
    }
    else
    {
        lastError = "Syntax Error: Third argument should be 'id' or 'var'!";
    }


    if(foundIndex)
    {
        if((atomIndices[0] < state_->atoms_.size()) && (atomIndices[1] < state_->atoms_.size())
           && (atomIndices[2] < state_->atoms_.size()))
        {
            *atom1Ptr = state_->atoms_[atomIndices[0]].get();
            *atom2Ptr = state_->atoms_[atomIndices[1]].get();
            *atom3Ptr = state_->atoms_[atomIndices[2]].get();

            succeeded = true;
        }
        else
        {
            lastError = "Invalid atom indices!";
        }
    }

    return std::pair<bool, std::string>(succeeded, lastError);
}

std::pair<bool, std::string> CParsingTools::retrieveBondAtoms(std::vector<std::string> arguments, size_t& arg, CAtom** atom1Ptr, CAtom** atom2Ptr, int* atomIndices) const
{
    std::string lastError;
    std::string text;
    bool foundIndex = false;
    bool succeeded = false;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "id")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[0] = atoi(text.data());
        text = CASCIIUtility::getArg(arguments, arg++);
        atomIndices[1] = atoi(text.data());
        foundIndex = true;
    }
    else if(text == "var")
    {
        int varIndex;

        text = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeBond))
        {
            CVarBond* p = (CVarBond*)varPtr;

            atomIndices[0] = p->atomIndex1_;
            atomIndices[1] = p->atomIndex2_;
            foundIndex = true;
        }
        else
        {
            lastError = "Variable missing or not of type 'bond'!";
        }
    }
    else
    {
        lastError = "Syntax Error: Third argument should be 'id' or 'var'!";
    }


    if(foundIndex)
    {
        if((atomIndices[0] < state_->atoms_.size()) && (atomIndices[1] < state_->atoms_.size()))
        {
            *atom1Ptr = state_->atoms_[atomIndices[0]].get();
            *atom2Ptr = state_->atoms_[atomIndices[1]].get();

            succeeded = true;
        }
        else
        {
            lastError = "Invalid atom indices!";
        }
    }

    return std::pair<bool, std::string>(succeeded, lastError);
}
