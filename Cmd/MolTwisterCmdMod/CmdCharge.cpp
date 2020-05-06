#include "CmdCharge.h"
#include "../../Utilities/ASCIIUtility.h"

namespace MolTwisterCmdMod
{
    std::string CCmdCharge::getCmd()
    {
        return "charge";
    }

    std::vector<std::string> CCmdCharge::getCmdLineKeywords()
    {
        return { "charge", "id", "var", "name" "to" };
    }

    std::vector<std::string> CCmdCharge::getCmdHelpLines()
    {
        return {
                    "charge id <atom index> to <partial charge>",
                    "charge var <variable name> to <partial charge>",
                    "charge name <atom ID> to <partial charge>"
               };
    }

    std::string CCmdCharge::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tModifies the charge of atoms given by\r\n";
        text+= "\t* the atom index, <atom index>\r\n";
        text+= "\t* the atom index contained in the variable <variable name>\r\n";
        text+= "\t* all atoms with atom ID = <atom ID> (e.g., O, H, C7)\r\n";
        text+= "\tto the charge given by <partial charge>.";

        return text;
    }

    std::string CCmdCharge::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text1, text2;
        CAtom* atomPtr = nullptr;
        bool bFoundIndex = false;
        std::vector<int> atomIndices;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 == "id")
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndices.emplace_back(atoi(text1.data()));
            bFoundIndex = true;
        }
        else if(text1 == "var")
        {
            int varIndex;

            text1 = CASCIIUtility::getArg(arguments, arg++);
            CVar* varPtr = state_->getVariable(text1.data(), varIndex);
            if(varPtr && (varPtr->getType() == CVar::typeAtom))
            {
                CVarAtom* p = (CVarAtom*)varPtr;

                atomIndices.emplace_back(p->atomIndex_);
                bFoundIndex = true;
            }
            else
            {
                lastError_ = "Variable missing or not of type 'atom'!";
                return lastError_;
            }
        }
        else if(text1 == "name")
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);

            for(int i=0; i<state_->atoms_.size(); i++)
            {
                text2 = state_->atoms_[i]->getID();
                if(text2 == text1)
                {
                    atomIndices.emplace_back(i);
                    bFoundIndex = true;
                }
            }
        }
        else
        {
            lastError_ = "Syntax Error: Third argument should be 'id', 'var' or 'name'!";
            return lastError_;
        }

        if(bFoundIndex)
        {
            double Q;

            text1 = CASCIIUtility::getArg(arguments, arg++);
            if(text1 != "to")
            {
                lastError_ = "Syntax Error: fifth argument should be 'to'!";
                return lastError_;
            }

            text1 = CASCIIUtility::getArg(arguments, arg++);
            Q = atof(text1.data());

            for(int i=0; i<atomIndices.size(); i++)
            {
                int atomIndex = atomIndices[i];

                if(atomIndex < state_->atoms_.size())
                {
                    atomPtr = state_->atoms_[atomIndex].get();
                    if(atomPtr)
                    {
                        atomPtr->Q_ = Q;
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
        }

        return lastError_;
    }
}
