#include "CmdRun.h"
#include "../../Utilities/ASCIIUtility.h"
#include "Simulator/MDSimulator.h"

#if INCLUDE_CUDA_COMMANDS == 1
namespace mtdev
{
    #include "Simulator/MDSimulator.h"
}
#endif

std::string CCmdRun::getCmd()
{
    return "run";
}

std::vector<std::string> CCmdRun::getCmdLineKeywords()
{
    return { "run", "cpu" };
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

    // Determine if simulation should be executed on CPU or on GPU
    std::string text = CASCIIUtility::getArg(arguments, arg++);
    bool runOnGPU = false;
    bool compiledForGPU = false;
    #if INCLUDE_CUDA_COMMANDS == 1
    {
        compiledForGPU = true;
    }
    #endif

    if(text == "cpu" )
    {
        runOnGPU = false;
        if(compiledForGPU)
        {
            fprintf(stdOut_, "\tMolTwister is compiled for GPU, but 'cpu' was explicitly specified. Hence, simulation will run on CPU!\r\n");
        }
        else
        {
            fprintf(stdOut_, "\tMolTwister is NOT compiled for GPU, but 'cpu' was explicitly specified. Hence, simulation will run on CPU!\r\n");
        }
    }
    else
    {
        if(compiledForGPU)
        {
            fprintf(stdOut_, "\tMolTwister is compiled for GPU and will run simulation on GPU!\r\n");
            runOnGPU = true;
        }
        else
        {
            fprintf(stdOut_, "\tMolTwister is NOT compiled for GPU and will therefore run simulation on CPU!\r\n");
            runOnGPU = false;
        }
    }

    // Run simulation
    if(runOnGPU)
    {
        #if INCLUDE_CUDA_COMMANDS == 1
        {
            std::stringstream stateContent;
            state_->serialize(stateContent, true);
            mtdev::CMDSimulator::run(1000000, 400, 4, "out.txt", stdOut_, stateContent, CMDSimulator::ensambleNPT);
        }
        #else
        {
            lastError_ = "\tMolTwister is NOT compiled for GPU, but the simulation attemptet to run on GPU!\r\n";
            return lastError_;
        }
        #endif
    }
    else
    {
        CMDSimulator::run(1000000, 400, 4, "out.txt", stdOut_, (CMolTwisterState*)state_, CMDSimulator::ensambleNPT);
    }

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

