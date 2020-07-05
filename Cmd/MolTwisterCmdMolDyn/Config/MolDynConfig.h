#pragma once
#include <stdio.h>
#include <vector>
#include <string>
#include "MolDynConfigStruct.h"

class CMolDynConfig
{
public:
    CMolDynConfig();

public:
    void print(FILE* stdOut=stdout);
    std::vector<std::string> getKeyWords();
    std::string set(std::string parameter, std::string value);
    void resetToDefaults();

public:
    SMolDynConfigStruct cfg_;
};
