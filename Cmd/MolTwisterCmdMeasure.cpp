#include "MolTwisterCmdMeasure.h"

#include "MolTwisterCmdMeasure/CmdBondLength.h"
#include "MolTwisterCmdMeasure/CmdBondLengthDyn.h"
#include "MolTwisterCmdMeasure/CmdBondSep.h"
#include "MolTwisterCmdMeasure/CmdAngle.h"
#include "MolTwisterCmdMeasure/CmdDihedral.h"
#include "MolTwisterCmdMeasure/CmdCoulombEnergy.h"
#include "MolTwisterCmdMeasure/CmdCoulombPotential.h"
#include "MolTwisterCmdMeasure/CmdAtomPos.h"
#include "MolTwisterCmdMeasure/CmdCenter.h"
#include "MolTwisterCmdMeasure/CmdPBC.h"
#include "MolTwisterCmdMeasure/CmdCharge.h"
#include "MolTwisterCmdMeasure/CmdCount.h"
#include "MolTwisterCmdMeasure/CmdBondCount.h"
#include "MolTwisterCmdMeasure/CmdOverlap.h"

using namespace MolTwisterCmdMeasure;

void CCmdMeasure::onAddKeywords()
{
    addKeyword("measure");
    addKeywords(parser_->getCmdLineKeywords());
}

void CCmdMeasure::onRegisterSubCommands()
{
    parser_->purge();
    parser_->registerCmd(std::make_shared<CCmdAngle>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdAtomPos>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdBondCount>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdBondLength>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdBondLengthDyn>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdBondSep>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCenter>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCharge>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCoulombEnergy>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCoulombPotential>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdCount>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdDihedral>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdOverlap>(state_, stdOut_));
    parser_->registerCmd(std::make_shared<CCmdPBC>(state_, stdOut_));
}

std::string CCmdMeasure::getHelpString() const
{
    return parser_->genHelpText(getCmd());
}

std::string CCmdMeasure::getHelpString(std::string subCommand) const
{
    return parser_->genHelpText(getCmd(), subCommand);
}

void CCmdMeasure::execute(std::string commandLine)
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
