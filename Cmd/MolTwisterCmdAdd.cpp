#include <iostream>
#include <math.h>
#include <vector>
#include "MolTwisterCmdAdd.h"
#include "MolTwisterCmdAdd/CmdAtom.h"
#include "MolTwisterCmdAdd/CmdAtoms.h"
#include "MolTwisterCmdAdd/CmdMDNonBonded.h"
#include "MolTwisterCmdAdd/CmdMDBond.h"
#include "MolTwisterCmdAdd/CmdMDAngle.h"
#include "MolTwisterCmdAdd/CmdMDDihedral.h"

void CCmdAdd::onAddKeywords()
{
    addKeyword("add");

    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdAdd::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdAtom>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdAtoms>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMDNonBonded>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMDBond>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMDAngle>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdMDDihedral>(state_, stdOut_));
}

std::string CCmdAdd::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdAdd::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdAdd::execute(std::string commandLine)
{
    std::string lastError = parser_->executeCmd(commandLine);
    if(!lastError.empty())
    {
        printf("%s\r\n", lastError.data());
        return;
    }
    
    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1)  state_->view3D_->requestUpdate(true);
        else                            state_->view3D_->requestUpdate(false);
    }    
}
