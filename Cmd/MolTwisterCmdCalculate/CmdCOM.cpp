#include "CmdCOM.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdCOM::getCmd()
{
    return "com";
}

std::vector<std::string> CCmdCOM::getCmdLineKeywords()
{
    return { "com", "sel" };
}

std::vector<std::string> CCmdCOM::getCmdHelpLines()
{
    return {
                "com <selection>"
           };
}

std::string CCmdCOM::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the center of mass (COM) of the specified selection. The possible selections are:\r\n";
    text+= "\t* sel - calculate the COM of the selected atoms.";

    return text;
}

std::string CCmdCOM::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string type;
    std::vector<int> atomIndices;

    type = CASCIIUtility::getArg(arguments, arg++);

    if(type == "sel")
    {
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            if(state_->atoms_[i]->isSelected()) atomIndices.emplace_back(i);
        }
    }
    else
    {
        lastError_ = "Error: expected 'sel'!";
        return lastError_;
    }

    C3DVector Rc = CMolTwisterStateTools(state_, stdOut_).getCenterOfMass(atomIndices, state_->currentFrame_);

    printf("\r\n");
    fprintf(stdOut_, "\tCenter of mass = (%.6f, %.6f, %.6f)\r\n", Rc.x_, Rc.y_, Rc.z_);

    return lastError_;
}
