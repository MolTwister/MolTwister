#pragma once
#include <string>

class CFileUtility
{
public:
    CFileUtility() {}
    
public:
    static void removeTabsFromFile(std::string fileName);
    static void addMissingSlash(std::string& folder);
    static std::string getCWD();
    static void removeExtension(std::string& filePath);
    static std::string getExtension(std::string filePath);
    static std::string readLine(FILE* fileHandle, bool& lastLineInFile);
    static std::string readToNextDelimiterIgnoreCppComments(FILE* fileHandle, char delimiter, bool& lastLineInFile);
    static std::string extractLeafDir(std::string absDir);
};
