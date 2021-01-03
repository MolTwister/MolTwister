#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "ASCIIUtility.h"

std::string CASCIIUtility::extractString(int startIndex, int endIndex, std::string stringIn)
{
    std::string stringOut;
    for(int i=startIndex; i<=endIndex; i++) stringOut+= stringIn[i];
    return stringOut;
}

bool CASCIIUtility::isWhiteSpace(char character, const char* whiteSpaceChars)
{
    int i=0;
    
    if(whiteSpaceChars != nullptr)
    {
        do
        {
            if(character == whiteSpaceChars[i++]) return true;
            
        } while(whiteSpaceChars[i] != '\0');
    }
    else
    {
        char szDefaultWhiteSpaces[3] = " \t";

        do
        {
            if(character == szDefaultWhiteSpaces[i++]) return true;
            
        } while(szDefaultWhiteSpaces[i] != '\0');
    }
    
    return false;
}

void CASCIIUtility::removeWhiteSpace(std::string& string, const char* whiteSpaceChars)
{
    std::string newStr;
    
    for(int i=0; i<string.length(); i++)
    {
        if(isWhiteSpace(string[i], whiteSpaceChars)) continue;
        
        newStr+= string[i];
    }
    
    string = newStr;
}

std::string CASCIIUtility::trimFromLeft(const std::string& str, const char* whiteSpaceChars, int maxTrim)
{
    std::string newStr;

    int numTrimmedChars = 0;
    int startFrom = 0;
    for(int i=0; i<str.length(); i++)
    {
        if(isWhiteSpace(str[i], whiteSpaceChars))
        {
            numTrimmedChars++;
            if(numTrimmedChars <= maxTrim)
                continue;
        }
        startFrom = i;
        break;
    }

    for(int i=startFrom; i<str.length(); i++)
    {
        newStr+= str[i];
    }

    return newStr;
}

std::string CASCIIUtility::trimFromRight(const std::string& str, const char* whiteSpaceChars, int maxTrim)
{
    std::string newStr;

    int numTrimmedChars = 0;
    int endIndex = int(str.length()-1);
    for(int i=int(str.length()-1); i>=0; i--)
    {
        if(isWhiteSpace(str[i], whiteSpaceChars))
        {
            numTrimmedChars++;
            if(numTrimmedChars <= maxTrim)
                continue;
        }
        endIndex = i;
        break;
    }

    for(int i=0; i<=endIndex; i++)
    {
        newStr+= str[i];
    }

    return newStr;
}

std::string CASCIIUtility::getWord(std::string line, int wordIndex, const char* whiteSpaceChars)
{
    size_t lineLen = line.size();
    bool inAWord = false;
    int currIndex;

    if(isWhiteSpace(line[0], whiteSpaceChars))
    {
        inAWord = false;
        currIndex = -1;
    }
    else
    {
        inAWord = true;
        currIndex = 0;
    }

    std::string word;
    for(int i=0; i<lineLen; i++)
    {
        if(inAWord)
        {
            if(isWhiteSpace(line[i], whiteSpaceChars))
            {
                inAWord = false;
            }
            else
            {
                if(wordIndex == currIndex)
                {
                    word+= line[i];
                }
            }
        }
        else
        {
            if(!isWhiteSpace(line[i], whiteSpaceChars))
            {
                inAWord = true;
                currIndex++;
                
                if(wordIndex == currIndex)
                {
                    word+= line[i];
                }
            }
        }
    }

    return word;
}

std::vector<std::string> CASCIIUtility::getWords(std::string line, const char* whiteSpaceChars)
{
    int wordIndex = 0;
    bool foundWord;
    std::string word;
    std::vector<std::string> words;

    do
    {
        word = getWord(line, wordIndex++, whiteSpaceChars);
        foundWord = (word.size() > 0) ? true : false;
        if(foundWord) words.emplace_back(word);

    } while(foundWord && (wordIndex < 10000000));

    return words;
}

std::vector<std::string> CASCIIUtility::getLines(std::string text)
{
    removeWhiteSpace(text, "\r");

    std::vector<std::string> lines;
    int len = (int)text.size();
    std::string line;
    for(int i=0; i<len; i++)
    {
        if(text[i] == '\n')
        {
            lines.emplace_back(line);
            line.clear();
        }
        else
        {
            line+= text[i];
        }
    }

    if(!line.empty()) lines.emplace_back(line);

    return lines;
}

std::string CASCIIUtility::getDelimitedWord(std::string line, int wordIndex, char startDelimiter, char endDelimiter)
{
    std::string word;
    size_t lineLen = line.size();
    bool inAWord = false;
    int currIndex = 0;
    int delimiterDepth = 0;

    for(int i=0; i<lineLen; i++)
    {
        if(line[i] == startDelimiter)
        {
            if(delimiterDepth == 0) inAWord = true;
            else
            {
                if(inAWord && (wordIndex == currIndex)) word+= line[i];
            }
            delimiterDepth++;
        }
        else if(line[i] == endDelimiter)
        {
            if(delimiterDepth == 1)
            {
                inAWord = false;
                currIndex++;
            }
            else
            {
                if(inAWord && (wordIndex == currIndex)) word+= line[i];
            }
            delimiterDepth--;
        }
        else
        {
            if(inAWord && (wordIndex == currIndex)) word+= line[i];
        }
        
        if(currIndex > wordIndex) return word;
    }
    
    return word;
}

std::string CASCIIUtility::removeCRLF(std::string text)
{
    std::string retString;
    size_t len = text.size();
    for(size_t i=0; i<len; i++)
    {
        if((text[i] != '\r') && (text[i] != '\n')) retString+= text[i];
    }

    return retString;
}

int CASCIIUtility::findString(std::string string, std::string line)
{
    size_t strLen = string.size();
    size_t lineLen = line.size();
    bool bMatch;
    
    for(int i=0; i<lineLen; i++)
    {
        bMatch = true;
        for(int j=0; j<strLen; j++)
        {
            if(((i + j) < lineLen) && (string[j] != line[i + j]))
                bMatch = false;
            if((i + j) >= lineLen)
                bMatch = false;
        }
        
        if(bMatch) return i;
    }
    
    return -1;
}

void CASCIIUtility::parseCLikeFuncCall(std::string inputString, std::string& funcName, std::string& arguments)
{
    int strLen = (int)inputString.size();
    int i;

    funcName.clear();
    arguments.clear();
    
    // Read function name
    for(i=0; i<strLen; i++)
    {
        if(inputString[i] == ' ') continue;
        if(inputString[i] == '(') break;
        
        funcName+= inputString[i];
    }
    
    // Read arguments
    int parenthesisLevel = 0;
    for(i=i+1; i<strLen; i++)
    {
        if((parenthesisLevel == 0) && (inputString[i] == ')')) break;
        
        if(inputString[i] == '(') parenthesisLevel++;
        if(inputString[i] == ')') parenthesisLevel--;
        
        arguments+= inputString[i];
    }
}

bool CASCIIUtility::isLineEmpty(std::string line)
{
    bool whiteSpace;
    char whiteSpaceDefs[] = " \t\r\n";
    int lineLen = (int)line.size();
    int numWhiteSpace = (int)strlen(whiteSpaceDefs);
    
    for(int i=0; i<lineLen; i++)
    {
        whiteSpace = false;
        for(int j=0; j<numWhiteSpace; j++)
            { if(line[i] == whiteSpaceDefs[j]) whiteSpace = true; }
        
        if(!whiteSpace) return false;
    }
    
    return true;
}

std::string CASCIIUtility::getArg(const std::vector<std::string>& arguments, size_t argIndex)
{
    if(argIndex < arguments.size())
    {
        return arguments[argIndex];
    }

    return "";
}

std::string CASCIIUtility::argsToString(const std::vector<std::string>& arguments)
{
    std::string argsString;

    bool first = true;
    for(auto argument : arguments)
    {
        if(!first) argsString+= " ";
        argsString+= argument;
        first = false;
    }

    return argsString;
}

std::string CASCIIUtility::argsToString(const std::vector<std::string>& arguments, size_t firstArgToInclude)
{
    std::vector<std::string> argumentsCpy = arguments;
    if(firstArgToInclude < argumentsCpy.size())
    {
        argumentsCpy.erase(argumentsCpy.begin(), argumentsCpy.begin() + firstArgToInclude);
        return CASCIIUtility::argsToString(argumentsCpy);
    }

    return "";
}

std::string CASCIIUtility::createMarkDownCodeBlock(std::string str, int numSpaces, bool removeFirstTab)
{
    std::vector<std::string> lines = getLines(str);

    std::string reconstructedString;
    reconstructedString+= "````text\r\n";
    for(std::string line : lines)
    {
        if(removeFirstTab) line = trimFromLeft(line, "\t", 1);
        line = trimFromRight(line);

        reconstructedString+= line;
        reconstructedString+= "\r\n";
    }
    reconstructedString+= "````\r\n";

    return reconstructedString;
}
