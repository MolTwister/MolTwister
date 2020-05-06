#include <iostream>
#include <math.h>
#include <float.h>
#include "Utilities/BashColor.h"
#include "MolTwisterCmd.h"

std::vector<std::string> CCmd::keywords_ = std::vector<std::string>();

void CCmd::addKeyword(std::string keyword)
{
    for(std::string cmpKey : keywords_)
    {
        if(cmpKey == keyword) return;
    }

    keywords_.emplace_back(keyword);
}

void CCmd::addKeywords(std::vector<std::string> keywords)
{
    for(auto keyword : keywords)
    {
        addKeyword(keyword);
    }
}

void CCmd::commonConstruct(CMolTwisterState* state)
{
    state_ = state;
    stdOut_ = stdout;
    skipTabRemoveReqFlag_ = false;
}
