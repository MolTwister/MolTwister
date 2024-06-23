#include <pthread.h>
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <Python.h>
#include <stdio.h>
#include "MolTwisterRestAPI.h"

CMolTwisterRestAPI::CMolTwisterRestAPI(const std::vector<std::string>& possibleFolderLocations)
{
    pythonMainFile_ = findPythonMainFile(possibleFolderLocations);
}

void CMolTwisterRestAPI::run() const
{
    pthread_t threadHandle;

    if(pthread_create(&threadHandle, nullptr, threadRun, (void*)this))
    {
        printf("Error: Could not create REST API thread!\r\n");
    }
}

void* CMolTwisterRestAPI::threadRun(void* arg)
{
    CMolTwisterRestAPI* ptrThis = (CMolTwisterRestAPI*)arg;
    std::string pythonString = readPythonFile(ptrThis->pythonMainFile_);

    PyGILState_STATE state = PyGILState_Ensure();
    PyRun_SimpleString(pythonString.data());
    PyGILState_Release(state);
    return nullptr;
}

std::string CMolTwisterRestAPI::findPythonMainFile(const std::vector<std::string>& possibleFolderLocations)
{
    std::string foundPythonFolder = "./Python/";

    for(const std::string& pythonFolder : possibleFolderLocations)
    {
        struct stat info;
        if((stat(pythonFolder.data(), &info) == 0) && (info.st_mode & S_IFDIR))
        {
            foundPythonFolder = pythonFolder;
            break;
        }
    }

    if((foundPythonFolder.size()) > 0 && (foundPythonFolder[foundPythonFolder.size()-1] != '/' ))
    {
        foundPythonFolder+= "/";
    }

    return foundPythonFolder + "MolTwisterAPI.py";
}

std::string CMolTwisterRestAPI::readPythonFile(const std::string& filePath)
{
    std::string document;
    std::ifstream file(filePath);

    if(file.is_open())
    {
        std::ostringstream stringStream;
        stringStream << file.rdbuf();
        document = stringStream.str();

        file.close();
    }

    return document;
}
