#include "CmdAtomPos.h"
#include "../../Utilities/ASCIIUtility.h"

namespace MolTwisterCmdMeasure
{
    std::string CCmdAtomPos::getCmd()
    {
        return "atompos";
    }

    std::vector<std::string> CCmdAtomPos::getCmdLineKeywords()
    {
        return { "atompos", "id", "var", "sel" };
    }

    std::vector<std::string> CCmdAtomPos::getCmdHelpLines()
    {
        return {
                    "atompos id <atom index>",
                    "atompos var <variable name>",
                    "atompos sel"
               };
    }

    std::string CCmdAtomPos::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the position of a single atom by\r\n";
        text+= "\t* defining the atom index directly by using the 'id' keyword\r\n";
        text+= "\t* obtaining atom indiex from a variable by using the 'var' keyword\r\n";
        text+= "\t* obtaining multiple atom indices from a selection by using the 'sel' keyword\r\n";
        text+= "\r\n";
        text+= "\tOutput (single atom measurement):\r\n";
        text+= "\tt(<x>, <y>, <z>)\r\n";
        text+= "\tOutput (multiple atom measurement):\r\n";
        text+= "\t1. [<atom index>, <atom ID>, <resname>]: (<x>, <y>, <z>)\r\n";
        text+= "\t2. [<atom index>, <atom ID>, <resname>]: (<x>, <y>, <z>)\r\n";
        text+= "\t     .\r\n";
        text+= "\t     .\r\n";
        text+= "\t     .\r\n";
        text+= "\tN. [<atom index>, <atom ID>, <resname>]: (<x>, <y>, <z>)\r\n";
        text+= "\twhere N is the number of selected atoms.";

        return text;
    }

    std::string CCmdAtomPos::execute(std::vector<std::string> arguments)
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
        else if(text1 == "sel")
        {
            fprintf(stdOut_, "\r\n");
            for(int i=0; i<state_->atoms_.size(); i++)
            {
                if(state_->atoms_[i]->isSelected())
                {
                    atomPtr = state_->atoms_[i].get();
                    text1 = atomPtr->getID();
                    if(!(state_->currentFrame_ < atomPtr->r_.size())) continue;
                    fprintf(stdOut_, "\t[%i,%s,%s]: (%.4f, %.4f, %.4f)\r\n", i, text1.data(), atomPtr->resname_.data(),
                            atomPtr->r_[state_->currentFrame_].x_, atomPtr->r_[state_->currentFrame_].y_, atomPtr->r_[state_->currentFrame_].z_);
                }
            }
            return lastError_;
        }
        else
        {
            lastError_ = "Syntax Error: Third argument should be 'id', 'var' or 'sel'!";
            return lastError_;
        }

        if(foundIndex)
        {
            if(atomIndex < state_->atoms_.size())
            {
                atomPtr = state_->atoms_[atomIndex].get();
                if(atomPtr)
                {
                    if(state_->currentFrame_ < atomPtr->r_.size())
                    {
                        fprintf(stdOut_, "\r\n\t(%.4f, %.4f, %.4f)\r\n", atomPtr->r_[state_->currentFrame_].x_, atomPtr->r_[state_->currentFrame_].y_, atomPtr->r_[state_->currentFrame_].z_);
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
}
