#include "MolTwisterCmdMolDyn.h"

#include "MolTwisterCmdMolDyn/CmdRun.h"
#include "MolTwisterCmdMolDyn/CmdFF.h"
#include "MolTwisterCmdMolDyn/CmdCfg.h"

#if INCLUDE_CUDA_COMMANDS == 1
#include "MolTwisterCmdMolDyn/CmdCudaTest.h"
#include "MolTwisterCmdMolDyn/CmdResetGPU.h"
#endif

void CCmdMolDyn::onAddKeywords()
{
    addKeyword("moldyn");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdMolDyn::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdRun>(state_, stdOut_, &molDynConfig_));
    parser_->registerCmd(std::make_shared<CCmdFF>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCfg>(state_, stdOut_, &molDynConfig_));

    #if INCLUDE_CUDA_COMMANDS == 1
    {
        parser_->registerCmd(std::make_shared<CCmdCudaTest>(state_, stdOut_));
        parser_->registerCmd(std::make_shared<CCmdResetGPU>(state_, stdOut_));
    }
    #endif
}

std::string CCmdMolDyn::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdMolDyn::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdMolDyn::execute(std::string commandLine)
{
    std::string lastError = parser_->executeCmd(commandLine);
    if(!lastError.empty())
    {
        printf("%s\r\n", lastError.data());
        return;
    }

    if(state_->view3D_)
    {
        if(state_->atoms_.size() == 1) state_->view3D_->requestUpdate(true);
        else                           state_->view3D_->requestUpdate(false);
    }
}
