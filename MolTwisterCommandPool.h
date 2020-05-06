#pragma once
#include "Cmd/MolTwisterCmd.h"

class CMolTwisterCommandPool
{
public:
    CMolTwisterCommandPool() {}
    
public:
    static void generateCmdList(CMolTwisterState* mtState, std::vector<std::shared_ptr<CCmd>>& cmdList);
};
