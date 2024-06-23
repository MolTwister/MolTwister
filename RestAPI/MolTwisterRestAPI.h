#pragma once
#include <string>
#include <vector>

class CMolTwisterRestAPI
{
public:
    CMolTwisterRestAPI(const std::vector<std::string>& possibleFolderLocations);

public:
    void run() const;

private:
    static void* threadRun(void* arg);
    static std::string findPythonMainFile(const std::vector<std::string>& possibleFolderLocations);
    static std::string readPythonFile(const std::string& filePath);

private:
    std::string pythonMainFile_;
};
