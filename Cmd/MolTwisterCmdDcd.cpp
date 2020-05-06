#include "MolTwisterCmdDcd.h"

#include "MolTwisterCmdDcd/CmdReadRecord.h"
#include "MolTwisterCmdDcd/CmdNumRecords.h"
#include "MolTwisterCmdDcd/CmdReadCoordinate.h"
#include "MolTwisterCmdDcd/CmdNumCoordinates.h"
#include "MolTwisterCmdDcd/CmdHeader.h"
#include "MolTwisterCmdDcd/CmdUnwrap.h"
#include "MolTwisterCmdDcd/CmdAtomicUnwrap.h"
#include "MolTwisterCmdDcd/CmdWrap.h"

void CCmdDcd::onAddKeywords()
{
    addKeyword("dcd");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdDcd::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdAtomicUnwrap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdHeader>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdNumCoordinates>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdNumRecords>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdReadCoordinate>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdReadRecord>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdUnwrap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdWrap>(state_, stdOut_));
}

std::string CCmdDcd::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdDcd::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdDcd::execute(std::string commandLine)
{
    std::string lastError = parser_->executeCmd(commandLine);
    if(!lastError.empty())
    {
        printf("%s\r\n", lastError.data());
        return;
    }
}
