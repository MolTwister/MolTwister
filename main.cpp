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
#include <stdio.h>
#include <stdlib.h>
#include "MolTwisterPythonExtensions.h"
#include "3DView/MolTwister3DView.h"
#include "RestAPI/MolTwisterRestAPI.h"

int main(int argc, char* argv[])
{
    CMolTwister mt;
    CMolTwisterRestAPI restAPI({ "./moltwister-api/Python", "/usr/local/bin/moltwister-api/Python" });
    C3DView view3D(argc, argv);
    
    // Initialize and configure Python interpreter
    // (also extend Python to include MolTwister function calls)
    g_pMT = &mt;

    #if PY_MAJOR_VERSION >= 3
    wchar_t name[] = L"moltwister";
    Py_SetProgramName(name);
    PyImport_AppendInittab("moltwister", &PyInitMolTwister);
    #else
    Py_SetProgramName(argv[0]);
    #endif
    Py_Initialize();
    PyEval_InitThreads();
    #if PY_MAJOR_VERSION < 3
    Py_InitModule("moltwister", MolTwisterMethods);
    #endif
    PyEval_SaveThread();

    // Set default host IP and port for the REST API
    std::string hostIP = "127.0.0.1";
    std::string hostPort = "5000";

    // Set headless mode and REST API to default be off
    bool headless = false;
    bool enableRestAPI = false;

    // Parse command line arguments
    std::vector<std::string> arguments;
    for(int i=0; i<argc; i++)
    {
        if(strcmp(argv[i], "--version") == 0)
        {
            printf("\r\n\tMolTwister version %s\r\n\r\n", MOLTWISTER_VER);
            Py_Finalize();
            return 0;
        }
        
        if(strcmp(argv[i], "-dirWorking") == 0) mt.startDir_ = CMolTwister::dirWorking;
        if(strcmp(argv[i], "-dirHome") == 0) mt.startDir_ = CMolTwister::dirHome;
        if(strcmp(argv[i], "-dirInitFile") == 0) mt.startDir_ = CMolTwister::dirInitFile;

        if(strcmp(argv[i], "-headless") == 0) headless = true;
        if(strcmp(argv[i], "-restapi") == 0) enableRestAPI = true;

        if(strcmp(argv[i], "--help") == 0)
        {
            printf("\r\n\tCommand line switches\r\n\t-------------------------------\r\n");
            printf("\t--version: display version number\r\n");
            printf("\t-dirWorking: open MolTwister in working directory (default)\r\n");
            printf("\t-dirHome: open MolTwister in home directory\r\n");
            printf("\t-dirInitFile: open MolTwister in directory specified in MolTwister.shortcuts.\r\n");
            printf("\t\tIf default entry does not exist, MolTwister is opened in home directory.\r\n");
            printf("\t-ip <host ip>: change host IP address for REST API. Default is 127.0.0.1.\r\n");
            printf("\t-port <host port>: change host port for REST API. Default is 5000.\r\n");
            printf("\t-headless: enable headless mode, where the 3D viewport is hidden.\r\n");
            printf("\t-restapi: enable REST API.\r\n");
            printf("\r\n");
            Py_Finalize();
            return 0;
        }

        if(strcmp(argv[i], "-ip") == 0)
        {
            i++;
            if(i < argc) hostIP = argv[i];
        }

        if(strcmp(argv[i], "-port") == 0)
        {
            i++;
            if(i < argc) hostPort = argv[i];
        }
    }

    // Start MolTwister command loop and 3D View
    mt.run(&view3D, hostIP, hostPort);
    if(enableRestAPI) restAPI.run();
    view3D.show(mt.getAtomsPtr(), mt.getGLObjsPtr(), mt.getCurrFramePtr(), mt.getDefaultAtProp(), headless);

    // Cleanup
    Py_Finalize();
    return 0;
}
