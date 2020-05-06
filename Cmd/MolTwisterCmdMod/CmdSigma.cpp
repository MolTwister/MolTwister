#include "CmdSigma.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolecularTools.h"

std::string CCmdSigma::getCmd()
{
    return "sigma";
}

std::vector<std::string> CCmdSigma::getCmdLineKeywords()
{
    return { "sigma", "id", "var", "to" };
}

std::vector<std::string> CCmdSigma::getCmdHelpLines()
{
    return {
                "sigma id <atom index> to <sigma>",
                "sigma var <variable name> to <sigma>"
           };
}

std::string CCmdSigma::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tChanges the van der Waals radius, denoted as sigma, to the value stated by <sigma>.\r\n";
    text+= "\tThis is not to be confued with the Lennard-Jones sigma for the assigned force fields.\r\n";
    text+= "\tThe van der Waals radius is used in calculations such as for iso surfaces describing\r\n";
    text+= "\ta surface plot surrounding a molecule or atomic cluster.\r\n";
    text+= "\r\n";
    text+= "\tThe atom index to receive the new value of sigma is identified by\r\n";
    text+= "\t* directly assigning the atom index, <atom index>\r\n";
    text+= "\t* a variable, <variable name>, containing the atom index";

    return text;
}

std::string CCmdSigma::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text1, text2;
    CAtom* atomPtr = nullptr;
    bool foundIndex = false;
    int atomIndex=0;

    text1 = CASCIIUtility::getArg(arguments, arg++);
    if(text1 == "id")
    {
        text1 = CASCIIUtility::getArg(arguments, arg++);
        atomIndex = atoi(text1.data());
        foundIndex = true;
    }
    else if(text1 == "var")
    {
        int varIndex;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        CVar* varPtr = state_->getVariable(text1.data(), varIndex);
        if(varPtr && (varPtr->getType() == CVar::typeAtom))
        {
            CVarAtom* p = (CVarAtom*)varPtr;

            atomIndex = p->atomIndex_;
            foundIndex = true;
        }
        else
        {
            lastError_ = "Variable missing or not of type 'atom'!";
            return lastError_;
        }
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should be 'id' or 'var'!";
        return lastError_;
    }

    if(foundIndex)
    {
        if(atomIndex < state_->atoms_.size())
        {
            atomPtr = state_->atoms_[atomIndex].get();
            if(atomPtr)
            {
                double sigma;

                text1 = CASCIIUtility::getArg(arguments, arg++);
                text2 = text1;

                text1 = CASCIIUtility::getArg(arguments, arg++);
                sigma = atof(text1.data());

                if(text2 == "to")
                {
                    atomPtr->sigma_ = sigma;
                }
                else
                {
                    lastError_ = "Syntax Error: fifth argument should be 'to'!";
                    return lastError_;
                }
            }
            else
            {
                lastError_ = "Could not find requested atom!";
                return lastError_;
            }
        }
        else
        {
            lastError_ = "Invalid atom index!";
            return lastError_;
        }
    }

    return lastError_;
}

