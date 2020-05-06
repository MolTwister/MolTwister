#pragma once
#include <string>
#include <vector>

class CASCIIUtility
{
public:
    static std::string extractString(int startIndex, int endIndex, std::string stringIn);
    static bool isWhiteSpace(char character, const char* whiteSpaceChars=nullptr);
    static void removeWhiteSpace(std::string& string, const char* whiteSpaceChars=nullptr);
    static std::string getWord(std::string line, int wordIndex, const char* whiteSpaceChars=nullptr);
    static std::vector<std::string> getWords(std::string line, const char* whiteSpaceChars=nullptr);
    static std::string getDelimitedWord(std::string line, int wordIndex, char startDelimiter, char endDelimiter);
    static std::string removeCRLF(std::string text);
    static int findString(std::string string, std::string line);
    static void parseCLikeFuncCall(std::string inputString, std::string& funcName, std::string& arguments);
    static bool isLineEmpty(std::string line);
    static std::string getArg(const std::vector<std::string>& arguments, size_t argIndex);
    static std::string argsToString(const std::vector<std::string>& arguments);
    static std::string argsToString(const std::vector<std::string>& arguments, size_t firstArgToInclude);
};
