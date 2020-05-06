#include "CmdAtomName.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/ConditionalOnXYZ.h"

std::string CCmdAtomName::getCmd()
{
    return "atomname";
}

std::vector<std::string> CCmdAtomName::getCmdLineKeywords()
{
    return { "atomname", "within" };
}

std::vector<std::string> CCmdAtomName::getCmdHelpLines()
{
    return {
                "atomname <atom ID to change> <atom ID to change to> [within <domain>]"
           };
}

std::string CCmdAtomName::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tRenames atom ID (e.g., H, O, C7) from <atom ID to change> to <atom ID to change to>, either\r\n";
    text+= "\tfor all atoms of the loaded system, if 'within' is not specified, or for the atoms within the\r\n";
    text+= "\tgiven <domain>.\r\n";
    text+= "\r\n";
    text+= CConditionalOnXYZ::getHelpText("<domain>", false);

    return text;
}

std::string CCmdAtomName::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CConditionalOnXYZ cond;
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
    }

    return lastError_;
}

