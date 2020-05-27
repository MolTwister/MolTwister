#include "CmdRun.h"
#include "../../Utilities/ASCIIUtility.h"
#include "MDLoop/MDLoop.h"
#include "MDLoop/Printf.h"

std::string CCmdRun::getCmd()
{
    return "run";
}

std::vector<std::string> CCmdRun::getCmdLineKeywords()
{
    return { "run" };
}

std::vector<std::string> CCmdRun::getCmdHelpLines()
{
    return {
        "run <atom ID to change> <atom ID to change to> [within <domain>]"
    };
}

std::string CCmdRun::getCmdFreetextHelp()
{
    std::string text;

    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";
    text+= "\t\r\n";

    return text;
}

std::string CCmdRun::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;

    int             iNStep = 50000;
    int             iOutputEvery = 400;
    int             iM = 4;
    std::string     szSystem;
    CMDLoop         MDLoop;
    CSimulationBox  SimBox(state_, stdOut_);

    SimBox.bNPTEnsemble = true;
    COut::SetOutputFile(fopen("out.txt", "w"));
    SimBox.InitSystem(CSimulationBox::sysLJCH4NormDens, iM);
    MDLoop.RunSimulation(SimBox, iNStep, iOutputEvery);
    COut::CloseOutputFile();

/*    CConditionalOnXYZ cond;
    std::string fromName, toName, text;
    bool within = false;

    fromName = CASCIIUtility::getArg(arguments, arg++);
    toName = CASCIIUtility::getArg(arguments, arg++);

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "within")
    {
        // Domain specifier of the form x<5&x>8|y<=4&y>3...
        // without any space in-between
        text = CASCIIUtility::getArg(arguments, arg++);
        cond.parse(text);
        within = true;
    }

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        std::string ID = state_->atoms_[i]->getID();

        if(ID == fromName)
        {
            if(!within)
            {
                state_->atoms_[i]->setID(toName.data());
            }
            else
            {
                if(state_->currentFrame_ >= state_->atoms_[i]->r_.size()) continue;
                C3DVector v = state_->atoms_[i]->r_[state_->currentFrame_];
                if(cond.isFulfilled(v.x_, v.y_, v.z_))
                {
                    state_->atoms_[i]->setID(toName.data());
                }
            }
        }
    }*/

    return lastError_;
}

