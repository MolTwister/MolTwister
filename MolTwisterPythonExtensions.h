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

#pragma once
#include "Utilities/ASCIIUtility.h"
#include "Utilities/DCDFile.h"
#include "MolTwisterCommandPool.h"
#include "Cmd/MolTwisterCmdPrint.h"
#include "MolTwister.h"
#include "Cmd/Tools/ProgressBar.h"

CMolTwister* g_pMT = nullptr;
CProgressBar g_progBar;

static PyObject* moltwister_mt_exec(PyObject* , PyObject* args)
{
    const char* argLinePtr;
    std::string stringArgLine, stringCommand, stringRet;
       
    if(!PyArg_ParseTuple(args, "s", &argLinePtr)) return nullptr;
    stringArgLine = argLinePtr;
    
    std::vector<std::shared_ptr<CCmd>> cmdList;
    CMolTwisterCommandPool::generateCmdList(g_pMT->getCurrentState(), cmdList);
    
    CASCIIUtility::removeWhiteSpace(stringArgLine, "\"");
    stringCommand = CASCIIUtility::getWord(stringArgLine, 0);

    int pipeSymbIndex = CASCIIUtility::findString("> ", stringArgLine);
    
    for(int i=0; i<(int)cmdList.size(); i++)
    {
        if(cmdList[i] && cmdList[i]->checkCmd(stringCommand.data()))
        {
            std::string fileName = ".moltwister_temp_file";
            if(pipeSymbIndex != -1)
            {
                fileName = stringArgLine.substr(pipeSymbIndex+1, std::string::npos);
                CASCIIUtility::removeWhiteSpace(fileName);
            }

            FILE* fileStdOut = fopen(fileName.data(), "w");
            if(fileStdOut) cmdList[i]->redirectOutput(fileStdOut);
            
            cmdList[i]->execute(stringArgLine.data());
            i = (int)cmdList.size();
            
            if(fileStdOut)
            {
                fclose(fileStdOut);
                fileStdOut = fopen(fileName.data(), "r");
                
                fseek(fileStdOut, 0L, SEEK_END);
                unsigned long size = ftell(fileStdOut);
                fseek(fileStdOut, 0L, SEEK_SET);
                
                std::string fileName;
                fileName.resize(size + 1);
                fread((char*)fileName.data(), sizeof(char), (size_t)size, fileStdOut);
                fileName[size] = '\0';
                stringRet = fileName;
                fclose(fileStdOut);
            }

            CFileUtility::removeTabsFromFile(fileName);
        }
    }

    return Py_BuildValue("s", stringRet.data());
}

static PyObject* moltwister_mt_get_num_atoms(PyObject*, PyObject* args)
{
    if(!PyArg_ParseTuple(args, ":get_num_atoms")) return nullptr;
    
    return Py_BuildValue("i", g_pMT->getCurrentState()->atoms_.size());
}

static PyObject* moltwister_mt_get_atom_pos(PyObject*, PyObject* args)
{
    int index, indexAxis;

    if(!PyArg_ParseTuple(args, "ii", &index, &indexAxis)) return nullptr;
    
    if(index < 0) return Py_BuildValue("d", 0.0);
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("d", 0.0);
    
    if(indexAxis == 0)
        return Py_BuildValue("d", g_pMT->getCurrentState()->atoms_[index]->r_[g_pMT->getCurrentState()->currentFrame_].x_);
    if(indexAxis == 1)
        return Py_BuildValue("d", g_pMT->getCurrentState()->atoms_[index]->r_[g_pMT->getCurrentState()->currentFrame_].y_);
    if(indexAxis == 2)
        return Py_BuildValue("d", g_pMT->getCurrentState()->atoms_[index]->r_[g_pMT->getCurrentState()->currentFrame_].z_);
    
    return Py_BuildValue("d", 0.0);
}

static PyObject* moltwister_mt_get_atom_type(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("s", "");
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("s", "");
    
    std::string ID = g_pMT->getCurrentState()->atoms_[index]->getID();
    
    return Py_BuildValue("s", ID.data());
}

static PyObject* moltwister_mt_get_atom_mass(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("d", 0.0);
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("d", 0.0);
    
    return Py_BuildValue("d", g_pMT->getCurrentState()->atoms_[index]->m_);
}

static PyObject* moltwister_mt_get_atom_charge(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("d", 0.0);
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("d", 0.0);
    
    return Py_BuildValue("d", g_pMT->getCurrentState()->atoms_[index]->Q_);
}

static PyObject* moltwister_mt_get_atom_resname(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("s", "");
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("s", "");
    
    return Py_BuildValue("s", g_pMT->getCurrentState()->atoms_[index]->resname_.data());
}

static PyObject* moltwister_mt_get_atom_molindex(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("i", 0);
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("i", 0);
    
    return Py_BuildValue("i", g_pMT->getCurrentState()->atoms_[index]->getMolIndex());
}

static PyObject* moltwister_mt_is_atom_sel(PyObject*, PyObject* args)
{
    int index;
    
    if(!PyArg_ParseTuple(args, "i", &index)) return nullptr;
    
    if(index < 0) return Py_BuildValue("i", -1);
    if(index >= (int)g_pMT->getCurrentState()->atoms_.size()) return Py_BuildValue("i", -1);
    
    return Py_BuildValue("i", g_pMT->getCurrentState()->atoms_[index]->isSelected() ? 1 : 0);
}

static PyObject* moltwister_mt_create_xyz_file(PyObject*, PyObject* args)
{
    const char* filePath;

    if(!PyArg_ParseTuple(args, "s", &filePath)) return nullptr;

    FILE* file = fopen(filePath, "w");
    if(file)
    {
        fclose(file);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_append_to_xyz_file(PyObject*, PyObject* args)
{
    const double toAU = 1.8897261245651;
    const char* filePath;
    double boxSizeX = 0.0;
    double boxSizeY = 0.0;
    double boxSizeZ = 0.0;
    PyObject* atomCoordinates = nullptr;
    bool convertToAU = false;

    if(!PyArg_ParseTuple(args, "sdddpO!", &filePath, &boxSizeX, &boxSizeY, &boxSizeZ, &convertToAU, &PyList_Type, &atomCoordinates)) return nullptr;
    Py_ssize_t numAtoms = PyList_Size(atomCoordinates);
    if(numAtoms < 0) return nullptr;

    std::vector<std::pair<std::string, C3DVector>> cppAtomCoordinates(numAtoms);
    for(Py_ssize_t i=0; i<numAtoms; i++)
    {
        PyObject* coordinate = PyList_GetItem(atomCoordinates, i);
        Py_ssize_t numEntries = PyList_Size(coordinate);
        if(numEntries != 4) continue;

        C3DVector cppCoordinate;
        std::string atomType;

        PyObject* entry = PyList_GetItem(coordinate, 0);
        PyObject* repr = PyObject_Repr(entry);
        PyObject* str = PyUnicode_AsEncodedString(repr, "utf-8", "~E~");
        atomType = PyBytes_AS_STRING(str);
        Py_XDECREF(repr);
        Py_XDECREF(str);

        entry = PyList_GetItem(coordinate, 1);
        cppCoordinate.x_ = PyFloat_AsDouble(entry);

        entry = PyList_GetItem(coordinate, 2);
        cppCoordinate.y_ = PyFloat_AsDouble(entry);

        entry = PyList_GetItem(coordinate, 3);
        cppCoordinate.z_ = PyFloat_AsDouble(entry);

        cppAtomCoordinates[i] = { atomType, cppCoordinate };
    }

    FILE* file = fopen(filePath, "a+");
    if(file)
    {
        int numAtoms = (int)cppAtomCoordinates.size();
        fprintf(file, "\t%i\r\n\tXYZ file, PBC = (%g, %g, %g)\r\n", numAtoms, boxSizeX, boxSizeY, boxSizeZ);
        for(int i=0; i<numAtoms; i++)
        {
            std::string ID = cppAtomCoordinates[i].first;
            C3DVector r = cppAtomCoordinates[i].second;
            double x = r.x_;
            double y = r.y_;
            double z = r.z_;
            if(convertToAU) { x*= toAU; y*= toAU; z*= toAU;}
            fprintf(file, "\t%-5s% -13.3f% -13.3f% -13.3f\r\n", ID.data(), x, y, z);
        }

        fprintf(file, "\r\n");
        fclose(file);
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_create_dcd_file(PyObject*, PyObject* args)
{
    const char* filePath;
    int numTimeSteps = 0;
    int stride = 1;
    double timeStep = 1.0;
    int numAtoms = 1;

    if(!PyArg_ParseTuple(args, "siidi", &filePath, &numTimeSteps, &stride, &timeStep, &numAtoms)) return nullptr;

    CDCDFile::createDCDFileIfNotExists(filePath, numTimeSteps, stride, timeStep, numAtoms);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_append_to_dcd_file(PyObject*, PyObject* args)
{
    const char* filePath;
    double boxSizeX = 0.0;
    double boxSizeY = 0.0;
    double boxSizeZ = 0.0;
    PyObject* atomCoordinates = nullptr;

    if(!PyArg_ParseTuple(args, "sdddO!", &filePath, &boxSizeX, &boxSizeY, &boxSizeZ, &PyList_Type, &atomCoordinates)) return nullptr;
    Py_ssize_t numAtoms = PyList_Size(atomCoordinates);
    if(numAtoms < 0) return nullptr;

    std::vector<C3DVector> cppAtomCoordinates(numAtoms);
    for(Py_ssize_t i=0; i<numAtoms; i++)
    {
        PyObject* coordinate = PyList_GetItem(atomCoordinates, i);
        Py_ssize_t numEntries = PyList_Size(coordinate);
        if(numEntries != 3) continue;

        C3DVector cppCoordinate;

        PyObject* entry = PyList_GetItem(coordinate, 0);
        cppCoordinate.x_ = PyFloat_AsDouble(entry);

        entry = PyList_GetItem(coordinate, 1);
        cppCoordinate.y_ = PyFloat_AsDouble(entry);

        entry = PyList_GetItem(coordinate, 2);
        cppCoordinate.z_ = PyFloat_AsDouble(entry);

        cppAtomCoordinates[i] = cppCoordinate;
    }

    CDCDFile dcdFile;

    std::function<std::tuple<double, double, double>(const int&)> getAtomPos = [&cppAtomCoordinates](const int& atomIndex)
    {
        C3DVector r = cppAtomCoordinates[atomIndex];
        return std::tuple<double, double, double>(r.x_, r.y_, r.z_);
    };

    CDCDFile::appendToDCDFile(filePath, (int)numAtoms, { boxSizeX, boxSizeY, boxSizeZ }, getAtomPos);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_begin_progress(PyObject*, PyObject* args)
{
    const char* argLinePtr;
    
    if(!PyArg_ParseTuple(args, "s", &argLinePtr)) return nullptr;
    g_progBar.beginProgress(argLinePtr);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_update_progress(PyObject*, PyObject* args)
{
    int step, totSteps;
    
    if(!PyArg_ParseTuple(args, "ii", &step, &totSteps)) return nullptr;
    g_progBar.updateProgress(step, totSteps);
    
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject* moltwister_mt_end_progress(PyObject*, PyObject* args)
{
    if(!PyArg_ParseTuple(args, ":end_progress")) return nullptr;
    g_progBar.endProgress();
    
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef MolTwisterMethods[] =
{
    // The mt_.. style commands are there for backward compatibility
    {"mt_exec", moltwister_mt_exec, METH_VARARGS, "Submit commandstring to MolTwister as a string. Result is returned as a string."},
    {"mt_get_num_atoms", moltwister_mt_get_num_atoms, METH_VARARGS, "Get number of atoms in 3D view."},
    {"mt_get_atom_pos", moltwister_mt_get_atom_pos, METH_VARARGS, "Get atom position."},
    {"mt_get_atom_type", moltwister_mt_get_atom_type, METH_VARARGS, "Get atom type."},
    {"mt_get_atom_mass", moltwister_mt_get_atom_mass, METH_VARARGS, "Get atom mass."},
    {"mt_get_atom_charge", moltwister_mt_get_atom_charge, METH_VARARGS, "Get atom charge."},
    {"mt_get_atom_resname", moltwister_mt_get_atom_resname, METH_VARARGS, "Get atom resname."},
    {"mt_get_atom_molindex", moltwister_mt_get_atom_molindex, METH_VARARGS, "Get atom molecular index."},
    {"mt_is_atom_sel", moltwister_mt_is_atom_sel, METH_VARARGS, "Check if atom is selected 1:yes, 0:no."},
    {"mt_create_xyz_file", moltwister_mt_create_xyz_file, METH_VARARGS, "Create XYZ file."},
    {"mt_append_to_xyz_file", moltwister_mt_append_to_xyz_file, METH_VARARGS, "Append to XYZ file from array of [atom type, x, y, z] lists."},
    {"mt_create_dcd_file", moltwister_mt_create_dcd_file, METH_VARARGS, "Create DCD file with basic header information."},
    {"mt_append_to_dcd_file", moltwister_mt_append_to_dcd_file, METH_VARARGS, "Append to DCD file from array of atom coordinates."},
    {"mt_begin_progress", moltwister_mt_begin_progress, METH_VARARGS, "Shows initial progress bar with given text."},
    {"mt_update_progress", moltwister_mt_update_progress, METH_VARARGS, "Updates progress bar according to given step information."},
    {"mt_end_progress", moltwister_mt_end_progress, METH_VARARGS, "Finishes progress bar, showing as 100 percent complete."},

    // Below is an improved naming convention, since in Python the name of the import will identify the commands as being for MolTwister
    {"exec", moltwister_mt_exec, METH_VARARGS, "Submit commandstring to MolTwister as a string. Result is returned as a string."},
    {"get_num_atoms", moltwister_mt_get_num_atoms, METH_VARARGS, "Get number of atoms in 3D view."},
    {"get_atom_pos", moltwister_mt_get_atom_pos, METH_VARARGS, "Get atom position."},
    {"get_atom_type", moltwister_mt_get_atom_type, METH_VARARGS, "Get atom type."},
    {"get_atom_mass", moltwister_mt_get_atom_mass, METH_VARARGS, "Get atom mass."},
    {"get_atom_charge", moltwister_mt_get_atom_charge, METH_VARARGS, "Get atom charge."},
    {"get_atom_resname", moltwister_mt_get_atom_resname, METH_VARARGS, "Get atom resname."},
    {"get_atom_molindex", moltwister_mt_get_atom_molindex, METH_VARARGS, "Get atom molecular index."},
    {"is_atom_sel", moltwister_mt_is_atom_sel, METH_VARARGS, "Check if atom is selected 1:yes, 0:no."},
    {"create_xyz_file", moltwister_mt_create_xyz_file, METH_VARARGS, "Create XYZ file."},
    {"append_to_xyz_file", moltwister_mt_append_to_xyz_file, METH_VARARGS, "Append to XYZ file from array of [atom type, x, y, z] lists."},
    {"create_dcd_file", moltwister_mt_create_dcd_file, METH_VARARGS, "Create DCD file with basic header information."},
    {"append_to_dcd_file", moltwister_mt_append_to_dcd_file, METH_VARARGS, "Append to DCD file from array of atom coordinates."},
    {"begin_progress", moltwister_mt_begin_progress, METH_VARARGS, "Shows initial progress bar with given text."},
    {"update_progress", moltwister_mt_update_progress, METH_VARARGS, "Updates progress bar according to given step information."},
    {"end_progress", moltwister_mt_end_progress, METH_VARARGS, "Finishes progress bar, showing as 100 percent complete."},

    {nullptr, nullptr, 0, nullptr}
};

#if PY_MAJOR_VERSION >= 3
PyModuleDef s_PyInterfaceModule =
{
    PyModuleDef_HEAD_INIT, "moltwister", nullptr, -1, MolTwisterMethods,
    nullptr, nullptr, nullptr, nullptr
};

PyObject* PyInitMolTwister()
{
    return PyModule_Create(&s_PyInterfaceModule);
}
#endif
