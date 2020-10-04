#pragma once
#include "Utilities/ASCIIUtility.h"
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
    if(!PyArg_ParseTuple(args, ":mt_get_num_atoms")) return nullptr;
    
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
    if(!PyArg_ParseTuple(args, ":mt_end_progress")) return nullptr;
    g_progBar.endProgress();
    
    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef MolTwisterMethods[] =
{
    {"mt_exec", moltwister_mt_exec, METH_VARARGS, "Submit commandstring to MolTwister as a string. Result is returned as a string."},
    {"mt_get_num_atoms", moltwister_mt_get_num_atoms, METH_VARARGS, "Get number of atoms in 3D view."},
    {"mt_get_atom_pos", moltwister_mt_get_atom_pos, METH_VARARGS, "Get atom position."},
    {"mt_get_atom_type", moltwister_mt_get_atom_type, METH_VARARGS, "Get atom type."},
    {"mt_get_atom_mass", moltwister_mt_get_atom_mass, METH_VARARGS, "Get atom mass."},
    {"mt_get_atom_charge", moltwister_mt_get_atom_charge, METH_VARARGS, "Get atom charge."},
    {"mt_get_atom_resname", moltwister_mt_get_atom_resname, METH_VARARGS, "Get atom resname."},
    {"mt_get_atom_molindex", moltwister_mt_get_atom_molindex, METH_VARARGS, "Get atom molecular index."},
    {"mt_is_atom_sel", moltwister_mt_is_atom_sel, METH_VARARGS, "Check if atom is selected 1:yes, 0:no."},
    {"mt_begin_progress", moltwister_mt_begin_progress, METH_VARARGS, "Shows initial progress bar with given text."},
    {"mt_update_progress", moltwister_mt_update_progress, METH_VARARGS, "Updates progress bar according to given step information."},
    {"mt_end_progress", moltwister_mt_end_progress, METH_VARARGS, "Finishes progress bar, showing as 100 percent complete."},
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
