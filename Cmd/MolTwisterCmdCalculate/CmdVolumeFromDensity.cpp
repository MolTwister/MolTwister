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

#include "CmdVolumeFromDensity.h"
#include "../../Utilities/ASCIIUtility.h"
#include <math.h>

#define NA 6.022E23 // Avogadros constant in mol^{-1}

std::string CCmdVolumeFromDensity::getCmd()
{
    return "volumefromdensity";
}

std::vector<std::string> CCmdVolumeFromDensity::getCmdLineKeywords()
{
    return { "volumefromdensity" };
}

std::vector<std::string> CCmdVolumeFromDensity::getCmdHelpLines()
{
    return {
        "volumefromdensity <target density in kg/m^3> <atomic masses in g/mol> <num. molecules>"
    };
}

std::string CCmdVolumeFromDensity::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tGiven a target density in kg/m^3, a list of atomic masses (comma separated with no space) for\r\n";
    text+= "\ta single atom in the system, as well as the number of such molecules in the system, the required\r\n";
    text+= "\tvolume in Angstrom^3 is calculated, in addition to the required side length required if the system\r\n";
    text+= "\twas cubic.";

    return text;
}

std::string CCmdVolumeFromDensity::execute(std::vector<std::string> arguments)
{
    lastError_ = "";
    size_t arg = 0;

    // Get target density (in kg/m^3)
    double targetDensity = atof(CASCIIUtility::getArg(arguments, arg++).data());

    // Get atomic masses of a single molecule
    std::string text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    std::vector<std::string> atomicMassList = CASCIIUtility::getWords(text, ",");

    // Get number of molecules
    int moleculeCount = atof(CASCIIUtility::getArg(arguments, arg++).data());

    // Calculate weight of single molecule (in g/mol)
    double M = 0.0;
    for(std::string m : atomicMassList)
    {
        M+= atof(m.data());
    }

    // Calculate the weight of all molecules
    M*= double(moleculeCount);

    // Convert to kg
    M/= (NA * 1000.0);

    // Calculate required volume (in m^3)
    double V = M / targetDensity;

    // Convert to AA^3
    V*= 1.0E30;

    printf("\r\n");
    for(int i=0; i<(int)atomicMassList.size(); i++)
    {
        fprintf(stdOut_, "\tMass of atom %i of single molecule = %.4f g/mol\r\n", i, atof(atomicMassList[i].data()));
    }
    fprintf(stdOut_, "\r\n\tTarget density = %.4f kg/m^3\r\n", targetDensity);
    fprintf(stdOut_, "\tNumber of molecules = %i\r\n", moleculeCount);
    fprintf(stdOut_, "\tDestination volume = %.4f Angstrom^3\r\n", V);
    fprintf(stdOut_, "\tDestination length for cubic system = %.4f Angstrom\r\n\r\n", pow(V, 1.0/3.0));

    return lastError_;
}
