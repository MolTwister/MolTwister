#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "MolTwisterPythonExtensions.h"
#include "3DView/MolTwister3DView.h"

int main(int argc, char* argv[])
{
    CMolTwister mt;
    C3DView view3D(argc, argv);
    
    // Initialize and configure Python interpreter
    // (also extend Python to include MolTwister function calls)
    g_pMT = &mt;

    #if PY_MAJOR_VERSION >= 3
    Py_SetProgramName(L"moltwister");
    PyImport_AppendInittab("moltwister", &PyInitMolTwister);
    #else
    Py_SetProgramName(argv[0]);
    #endif
    Py_Initialize();
    #if PY_MAJOR_VERSION < 3
    Py_InitModule("moltwister", MolTwisterMethods);
    #endif
    
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

        if(strcmp(argv[i], "--help") == 0)
        {
            printf("\r\n\tCommand line switches\r\n\t-------------------------------\r\n");
            printf("\t--version: display version number\r\n");
            printf("\t-dirWorking: open MolTwister in working directory (default)\r\n");
            printf("\t-dirHome: open MolTwister in home directory\r\n");
            printf("\t-dirInitFile: open MolTwister in directory specified in MolTwister.shortcuts.\r\n");
            printf("\t\tIf default entry does not exist MolTwister is opened in home directory.\r\n");
            printf("\r\n");
            Py_Finalize();
            return 0;
        }
    }

    // Start MolTwister command loop and 3D View
    mt.run(&view3D);
    view3D.show(mt.getAtomsPtr(), mt.getGLObjsPtr(), mt.getCurrFramePtr(), mt.getDefaultAtProp());

    // Cleanup
    Py_Finalize();
    return 0;
}
