#include "CmdAtomPos.h"
#include "../../Utilities/ASCIIUtility.h"

namespace MolTwisterCmdMod
{
    std::string CCmdAtomPos::getCmd()
    {
        return "atompos";
    }

    std::vector<std::string> CCmdAtomPos::getCmdLineKeywords()
    {
        return { "atompos", "id", "var", "all", "sel", "to", "by", "geomcent", "flip" };
    }

    std::vector<std::string> CCmdAtomPos::getCmdHelpLines()
    {
        return {
                    "atompos id <atom index> to <x> <y> <z>",
                    "atompos id <atom index> by <x> <y> <z>",
                    "atompos id <atom index> geomcent <x> <y> <z>",
                    "atompos id <atom index> flip <flip axis>",
                    "atompos var <variable name> to <x> <y> <z>",
                    "atompos var <variable name> by <x> <y> <z>",
                    "atompos var <variable name> geomcent <x> <y> <z>",
                    "atompos var <variable name> flip <flip axis>",
                    "atompos all to <x> <y> <z>",
                    "atompos all by <x> <y> <z>",
                    "atompos all geomcent <x> <y> <z>",
                    "atompos all flip <flip axis>",
                    "atompos sel to <x> <y> <z>",
                    "atompos sel by <x> <y> <z>",
                    "atompos sel geomcent <x> <y> <z>",
                    "atompos sel flip <flip axis>"
               };
    }

    std::string CCmdAtomPos::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tModifies the atomic positions of one or more atoms. The atoms to be moved can be specified by\r\n";
        text+= "\t* id <atom index>\r\n";
        text+= "\t* var <variable name>, where variable name contains atom index\r\n";
        text+= "\t* all, which signifies all loaded or created atoms\r\n";
        text+= "\t* sel, which signifies all selected atoms\r\n";
        text+= "\tThe relative motion of the translations can be one of the following\r\n";
        text+= "\t* to <x> <y> <z>, which will move all chosen atoms to (<x>, <y>, <z>)\r\n";
        text+= "\t* by <x> <y> <z>, which will translate all chosen atoms by <x>, <y> and <z> along the x-, y-\r\n";
        text+= "\t                  and z-axis, respectively\r\n";
        text+= "\t* geomcent <x> <y> <z>, which will translate all chosen atoms such that geometric center of\r\n";
        text+= "\t                        their bounding box is positioned at (<x>, <y>, <z>)\r\n";
        text+= "\t* flip <flip axis>, where <flip axis> is either x, y, or z, which will flip all the chosen\r\n";
        text+= "\t                    atoms are mirrored across the selected axis around the bounding box center\r\n";
        text+= "\t                    of the chosen atoms";

        return text;
    }

    std::string CCmdAtomPos::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text1;
        CAtom* atomPtr = nullptr;
        std::vector<int> atomIndices;

        text1 = CASCIIUtility::getArg(arguments, arg++);
        if(text1 == "id")
        {
            text1 = CASCIIUtility::getArg(arguments, arg++);
            atomIndices.emplace_back(atoi(text1.data()));
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
            }
            else
            {
                lastError_ = "Variable missing or not of type 'atom'!";
                return lastError_;
            }
        }
        else if(text1 == "all")
        {
            for(int i=0; i<state_->atoms_.size(); i++)
            {
                atomIndices.emplace_back(i);
            }
        }
        else if(text1 == "sel")
        {
            for(int i=0; i<state_->atoms_.size(); i++)
            {
                if(state_->atoms_[i]->isSelected())
                    atomIndices.emplace_back(i);
            }
        }
        else
        {
            lastError_ = "Syntax Error: Third argument should be 'id', 'var' or 'all'!";
            return lastError_;
        }

        if(atomIndices.size() > 0)
        {
            std::vector<C3DRect> boundBoxes;
            C3DVector v, v0;
            int method = 0;
            int flipAxis = -1;

            text1 = CASCIIUtility::getArg(arguments, arg++);

            if(text1 == "to")               method = 1;
            else if(text1 == "by")          method = 2;
            else if(text1 == "geomcent")    method = 3;
            else if(text1 == "flip")        method = 4;
            else
            {
                lastError_ = "Syntax Error: fifth argument should be 'to', 'by', 'flip' or 'geomcent'!";
                return lastError_;
            }

            for(int i=0; i<atomIndices.size(); i++)
            {
                if((atomIndices[i] >= state_->atoms_.size()) || (atomIndices[i] < 0))
                {
                    lastError_ = std::string("Invalid atom index ") + std::to_string(atomIndices[i]) + std::string("!");
                    return lastError_;
                }

                atomPtr = state_->atoms_[atomIndices[i]].get();
                if(!atomPtr)
                {
                    lastError_ = std::string("Could not find requested atom (index ") + std::to_string(atomIndices[i]) + std::string(")!");
                    return lastError_;
                }
            }

            if((method == 1) || (method == 2) || (method == 3))
            {
                text1 = CASCIIUtility::getArg(arguments, arg++);
                v.x_ = atof(text1.data());
                text1 = CASCIIUtility::getArg(arguments, arg++);
                v.y_ = atof(text1.data());
                text1 = CASCIIUtility::getArg(arguments, arg++);
                v.z_ = atof(text1.data());
            }

            if(method == 3)
            {
                if(state_->atoms_.size() <= 0)
                {
                    lastError_ = "Error, no atoms in system!";
                    return lastError_;
                }

                atomPtr = state_->atoms_[0].get();
                for(int j=0; j<atomPtr->r_.size(); j++)
                    boundBoxes.emplace_back(state_->calcBoundingBox(j, atomIndices));
            }

            if(method == 4)
            {
                text1 = CASCIIUtility::getArg(arguments, arg++);
                if(text1 == "x")        flipAxis = 0;
                else if(text1 == "y")   flipAxis = 1;
                else if(text1 == "z")   flipAxis = 2;
                else
                {
                    lastError_ = "Error, possible flip axis are 'x', 'y', or 'z'!";
                    return lastError_;
                }

                if(state_->atoms_.size() <= 0)
                {
                    lastError_ = "Error, no atoms in system!";
                    return lastError_;
                }

                atomPtr = state_->atoms_[0].get();
                for(int j=0; j<atomPtr->r_.size(); j++)
                    boundBoxes.emplace_back(state_->calcBoundingBox(j, atomIndices));
            }

            for(int i=0; i<atomIndices.size(); i++)
            {
                atomPtr = state_->atoms_[atomIndices[i]].get();
                for(int j=0; j<atomPtr->r_.size(); j++)
                {
                    switch(method)
                    {
                        case 1:
                            atomPtr->r_[j] = v;
                            break;
                        case 2:
                            atomPtr->r_[j]+= v;
                            break;
                        case 3:
                            if(j >= boundBoxes.size()) continue;
                            v0.x_ = v.x_ - (boundBoxes[j].rLow_.x_ + boundBoxes[j].rHigh_.x_) / 2.0;
                            v0.y_ = v.y_ - (boundBoxes[j].rLow_.y_ + boundBoxes[j].rHigh_.y_) / 2.0;
                            v0.z_ = v.z_ - (boundBoxes[j].rLow_.z_ + boundBoxes[j].rHigh_.z_) / 2.0;
                            atomPtr->r_[j]+= v0;
                            break;
                        case 4:
                            if(j >= boundBoxes.size()) continue;
                            if(flipAxis == 0)
                            {
                                double xc = (boundBoxes[j].rLow_.x_ + boundBoxes[j].rHigh_.x_) / 2.0;
                                atomPtr->r_[j].x_ = 2.0*xc - atomPtr->r_[j].x_;
                            }
                            if(flipAxis == 1)
                            {
                                double yc = (boundBoxes[j].rLow_.y_ + boundBoxes[j].rHigh_.y_) / 2.0;
                                atomPtr->r_[j].y_ = 2.0*yc - atomPtr->r_[j].y_;
                            }
                            if(flipAxis == 2)
                            {
                                double zc = (boundBoxes[j].rLow_.z_ + boundBoxes[j].rHigh_.z_) / 2.0;
                                atomPtr->r_[j].z_ = 2.0*zc - atomPtr->r_[j].z_;
                            }
                            break;
                    }
                }
            }
        }

        return lastError_;
    }
}
