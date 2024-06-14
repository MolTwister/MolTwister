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

#include "MolTwisterCmdPython.h"
#include "MolTwisterCmdDcd.h"

void CCmdPython::onAddKeywords()
{
    addKeyword("mtpython");
}

std::string CCmdPython::getHelpString() const
{
    std::string text;

    text+= "\tUsage: mtpython {<Python line of code>}\r\n";
    text+= "\r\n";
    text+= "\tAny Python line of code can be executed, one by one. For example, the sequence:\r\n";
    text+= "\r\n";
    text+= "\tmtpython {a = 5}\r\n";
    text+= "\tmtpython {b = 20}\r\n";
    text+= "\tmtpython {print(\"Output: %f\" % (a+b))}\r\n";
    text+= "\r\n";
    text+= "\twould produce the output 'Output: 25.000000'. Note that Python 'import' will load\r\n";
    text+= "\tthe specified library and will be available for the next use of 'mtpython'.\r\n";
    text+= "\t\r\n";
    text+= "\tBy importing the 'moltwister' library (e.g., import moltwister as mt), several\r\n";
    text+= "\tPython functions will be made available that can query / manipulate the state of\r\n";
    text+= "\tMolTwister. These are as follows\r\n";
    text+= "\r\n";
    text+= "\t\texec(<moltwister command>) : result as string\r\n";
    text+= "\t\tExecutes a moltwister command. For example, exec(\"list all\"), which\r\n";
    text+= "\t\twill return result of 'list all'.\r\n";
    text+= "\r\n";
    text+= "\t\tget_num_atoms() : result as int\r\n";
    text+= "\t\tReturns the number of atoms presently loaded or added.\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_pos(atomIndex:integer, axisIndex:integer) : result as integer\r\n";
    text+= "\t\tReturns the atom position of a given atom index and a given axis (0=x, 1=y, 2=z).\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_type(atomIndex:integer) : result as string\r\n";
    text+= "\t\tReturns the atom type of the atom at the given atom index.\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_mass(atomIndex:integer) : result as integer\r\n";
    text+= "\t\tReturns the assigned atomic mass of the atom at the given atom index.\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_charge(atomIndex:integer) : result as integer\r\n";
    text+= "\t\tReturns the assigned atomic charge of the atom at the given atom index.\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_resname(atomIndex:integer) : result as string\r\n";
    text+= "\t\tReturns the assigned resname of the atom at the given atom index.\r\n";
    text+= "\r\n";
    text+= "\t\tget_atom_molindex(atomIndex:integer) : result as integer\r\n";
    text+= "\t\tReturns the assigned molecular index of the atom at the given atom index.\r\n";
    text+= "\r\n";
    text+= "\t\tis_atom_sel(atomIndex:integer) : result as boolean\r\n";
    text+= "\t\tReturns true if the atom at the given atom index is selected, else it returns false.\r\n";
    text+= "\r\n";
    text+= "\t\tcreate_xyz_file(filePath:string) : no result\r\n";
    text+= "\t\tCreate an empty XYZ file.\r\n";
    text+= "\r\n";
    text+= "\t\tappend_to_xyz_file(filePath:string, boxSizeX:float, boxSizeY:float, boxSizeZ:float, convertToAU:bool, atomCoordinates:list) : no result\r\n";
    text+= "\t\tAppend list of [atomTypeString, x, y, z]-lists to XYZ file.\r\n";
    text+= "\r\n";
    text+= "\t\tcreate_dcd_file(filePath:string, numTimeSteps:int, stride:int, timeStep:float, numAtoms:int) : no result\r\n";
    text+= "\t\tCreate a DCD file with given header information.\r\n";
    text+= "\r\n";
    text+= "\t\tappend_to_dcd_file(filePath:string, boxSizeX:float, boxSizeY:float, boxSizeZ:float, atomCoordinates:list) : no result\r\n";
    text+= "\t\tAppend list of [x, y, z]-lists to DCD file.\r\n";
    text+= "\r\n";
    text+= "\t\tbegin_progress(progBarDescription:string) : no result\r\n";
    text+= "\t\tShows initial progress bar (in the command line shell) with given text.\r\n";
    text+= "\r\n";
    text+= "\t\tupdate_progress(step:integer, totalSteps:integer) : no result\r\n";
    text+= "\t\tUpdates progress bar according to given step information.\r\n";
    text+= "\r\n";
    text+= "\t\tend_progress() : no result\r\n";
    text+= "\t\tFinishes progress bar, showing as 100 percent complete.\r\n";
    text+= "\t\r\n";
    text+= "\tAn example sequence could be:\r\n";
    text+= "\r\n";
    text+= "\tmtpython {import moltwister as mt}\r\n";
    text+= "\tmtpython {print(\"Num atoms: %f\" % mt.get_num_atoms())}\r\n";
    text+= "\tmtpython {print(\"%s\" % mt.exec(\"list all\"))}\r\n";
    text+= "\r\n";
    text+= "\tNote that all the above examples are written in the language of Python version 3.0\r\n";
    text+= "\tand may not be compatible with earlier versions of Python.\r\n";
    text+= "\t\r\n";
    text+= "\tIt is also possible to write a python script and then load this script using\r\n";
    text+= "\tthe 'load python <name of script>' command, where only the <Python line of code>\r\n";
    text+= "\tparts of the mtpython commands are included within the script.";

    return text;
}

void CCmdPython::execute(std::string commandLine)
{
    if(!state_) return;
    
    std::string pythonCodeSnippet;
    pythonCodeSnippet = CASCIIUtility::getDelimitedWord(commandLine, 0, '{', '}');

    /////////////////////////////////////////////////////////
    // More information available on:
    // https://docs.python.org/2/extending/embedding.html
    /////////////////////////////////////////////////////////
    
    PyGILState_STATE gilState = PyGILState_Ensure();
    PyRun_SimpleString(pythonCodeSnippet.data());
    PyGILState_Release(gilState);
}
