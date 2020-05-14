#include <iostream>
#include <vector>
#include "Utilities/BashColor.h"
#include "Utilities/3DRect.h"
#include "MolTwisterCmdGet.h"
#include "Tools/MolTwisterStateTools.h"

#if INCLUDE_CUDA_COMMANDS == 1
#include "Cmd/Tools/CudaDeviceList.h"
#endif

void CCmdGet::onAddKeywords()
{
    addKeyword("get");
    addKeyword("atomtypes");
    addKeyword("mdinconsistency");
    addKeyword("bondinfo");
    addKeyword("userdefpbc");
    addKeyword("gpuinfo");
}

std::string CCmdGet::getHelpString() const
{
    std::string  text;
    
    text+= "\tUsage: get <args>\r\n";
    text+= "\r\n";
    text+= "\tRetrieve various properties of the system. <args> can be:\r\n";
    text+= "\r\n";
    text+= "\t       * atomtypes                 :   List all atom types present in system\r\n";
    text+= "\t       * mdinconsistency           :   Check if any inconsistencies in MD force-field can be found\r\n";
    text+= "\t       * bondinfo <ind 1> <ind 2>  :   Get a list of bonds connected to atom indices 1 and 2\r\n";
    text+= "\t       * userdefpbc                :   Get user defined PBCs\r\n";
    text+= "\t       * gpuinfo                   :   Check if CUDA is available and return info\r\n";
    
    return text;
}

void CCmdGet::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    
    if(text == "atomtypes")
    {
        parseAtomtypesCommand(commandLine, arg);
    }
    
    else if(text == "mdinconsistency")
    {
        parseMdinconsistencyCommand(commandLine, arg);
    }

    else if(text == "bondinfo")
    {
        parseBondinfoCommand(commandLine, arg);
    }

    else if(text == "userdefpbc")
    {
        parseUserdefpbcCommand(commandLine, arg);
    }

    else if(text == "gpuinfo")
    {
        parseGpuinfoCommand(commandLine, arg);
    }
    
    else
    {
        printf("Syntax Error: Second argument should specify what to get!");
    }
}

void CCmdGet::parseAtomtypesCommand(std::string, int&)
{
    std::vector<std::string> listOfAtomTypes;
    std::vector<std::string> listOfResnames;
    CAtom* atomPtr;
    
    state_->searchForAtomTypes(listOfAtomTypes, &listOfResnames);
    
    printf("\r\n");
    for(int i=0; i<listOfAtomTypes.size(); i++)
    {
        atomPtr = state_->getFirstOccurenceOf(listOfAtomTypes[i]);
        fprintf(stdOut_, "\t%i\t %s m=%.4f q=%.4f resname=%s\r\n", i+1, listOfAtomTypes[i].data(), atomPtr ? atomPtr->m_ : 0.0, atomPtr ? atomPtr->Q_ : 0.0, listOfResnames[i].data());
    }
}

void CCmdGet::parseMdinconsistencyCommand(std::string, int&)
{
    CMolTwisterStateTools(state_, stdOut_).reportConsistencyOfMDForceField();
}

void CCmdGet::parseBondinfoCommand(std::string commandLine, int& arg)
{
    std::string text;
    int indexAtom1, indexAtom2;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    indexAtom1 = atoi(text.data());

    text = CASCIIUtility::getWord(commandLine, arg++);
    indexAtom2 = atoi(text.data());
    
    if((indexAtom1 < 0) || (indexAtom1 >= state_->atoms_.size()))
    {
        printf("Error: bond %i could not be found!\r\n", indexAtom1);
        return;
    }

    if((indexAtom2 < 0) || (indexAtom2 >= state_->atoms_.size()))
    {
        printf("Error: bond %i could not be found!\r\n", indexAtom1);
        return;
    }
    
    if(state_->currentFrame_ < 0)
    {
        printf("Error: found no frames!\r\n");
        return;
    }
    
    CAtom* atom1Ptr = state_->atoms_[indexAtom1].get();
    CAtom* atom2Ptr = state_->atoms_[indexAtom2].get();
    
    fprintf(stdOut_, "\r\n\tAtom %i bonded to:\r\n", indexAtom1);
    for(int i=0; i<atom1Ptr->getNumBonds(); i++)
    {
        double dR = -1.0;
        CAtom* atomPtr = atom1Ptr->getBondDest(i);
        if(atomPtr == atom2Ptr) CBashColor::setSpecial(CBashColor::specBright);
        if(state_->currentFrame_ < atomPtr->r_.size()) dR = atomPtr->getDistanceTo(atom1Ptr, state_->currentFrame_);
        fprintf(stdOut_, "\t* Atom %i -> R=%.6f\r\n", state_->getAtomIndex(atomPtr), dR);
        CBashColor::setSpecial();
    }

    fprintf(stdOut_, "\r\n\tAtom %i bonded to:\r\n", indexAtom2);
    for(int i=0; i<atom2Ptr->getNumBonds(); i++)
    {
        double dR = -1.0;
        CAtom* atomPtr = atom2Ptr->getBondDest(i);
        if(atomPtr == atom1Ptr) CBashColor::setSpecial(CBashColor::specBright);
        if(state_->currentFrame_ < atomPtr->r_.size()) dR = atomPtr->getDistanceTo(atom2Ptr, state_->currentFrame_);
        fprintf(stdOut_, "\t* Atom %i -> R=%.6f\r\n", state_->getAtomIndex(atomPtr), dR);
        CBashColor::setSpecial();
    }
}

void CCmdGet::parseUserdefpbcCommand(std::string, int&)
{
    if(!state_->view3D_)
    {
        printf("Error: Could not find 3D View!\r\n");
        return;
    }
    
    C3DRect pbc = state_->view3D_->getUserPBC();
    
    fprintf(stdOut_, "\r\n\tUser defined PBC %s\r\n", state_->view3D_->isUserPBCEnabled() ? "enabled" : "disabled");
    fprintf(stdOut_, "\tUser_x = [%.4f, %.4f]\r\n", pbc.rLow_.x_, pbc.rHigh_.x_);
    fprintf(stdOut_, "\tUser_y = [%.4f, %.4f]\r\n", pbc.rLow_.y_, pbc.rHigh_.y_);
    fprintf(stdOut_, "\tUser_z = [%.4f, %.4f]\r\n", pbc.rLow_.z_, pbc.rHigh_.z_);
}

void CCmdGet::parseGpuinfoCommand(std::string, int&)
{
    #if INCLUDE_CUDA_COMMANDS == 1
    {
        CCudaDeviceList devList;
        devList.print(stdOut_);
    }
    #else
    {
        fprintf(stdOut_, "\r\n\tNo CUDA acceleration available. You need to have an Nvidia compatible GPU, the CUDA toolkit installed and set INCLUDE_CUDA_COMMANDS=1 in CMake!\r\n");
    }
    #endif
}

