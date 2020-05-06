#include <iostream>
#include "MolTwisterCmdVar.h"

void CCmdVar::onAddKeywords()
{
    addKeyword("var");
    addKeyword("atom");
    addKeyword("bond");
    addKeyword("angle");
}

std::string CCmdVar::getHelpString() const
{
    std::string text;
    
    text+= "\tUsage: var <type+value> <variable name>\r\n";
    text+= "\r\n";
    text+= "\tDefine a variable that, for example, can be used (in some situations)\r\n";
    text+= "\tinstead of atomic indices. <type+value> can be one of the following:\r\n";
    text+= "\r\n";
    text+= "\t  * atom <N>                     : Define atomic variable with index <N>\r\n";
    text+= "\t  * bond <N1> <N2>               : Define bond variable between index <N1>\r\n";
    text+= "\t                                   and index <N2>\r\n";
    text+= "\t  * angle <N1> <N2> <N3>         : Define angle variable between index <N1>,\r\n";
    text+= "\t                                   <N2> and <N3>\r\n";
    text+= "\t  * dihedral <N1> <N2> <N3> <N4> : Define dihedral variable between index\r\n";
    text+= "\t                                   <N1>, <N2>, <N3> and <N4>\r\n";
    text+= "\r\n";
    text+= "\tDefining variables can be useful for commands such as 'mod', 'measure',\r\n";
    text+= "\tand 'gauss9', to avoid having to remember atomic indices.";
    
    return text;
}

void CCmdVar::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "atom")
    {
        CVarAtom varAtom;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varAtom.atomIndex_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varAtom.setName(text.data());
        
        state_->addVariable(varAtom);
    }
    else if(text == "bond")
    {
        CVarBond varBond;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varBond.atomIndex1_ = atoi(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        varBond.atomIndex2_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varBond.setName(text.data());
        
        state_->addVariable(varBond);
    }
    else if(text == "angle")
    {
        CVarAngle varAngle;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varAngle.atomIndex1_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varAngle.atomIndex2_ = atoi(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        varAngle.atomIndex3_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varAngle.setName(text.data());
        
        state_->addVariable(varAngle);
    }
    else if(text == "dihedral")
    {
        CVarDihedral varDihedral;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varDihedral.atomIndex1_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varDihedral.atomIndex2_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varDihedral.atomIndex3_ = atoi(text.data());

        text = CASCIIUtility::getWord(commandLine, arg++);
        varDihedral.atomIndex4_ = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        varDihedral.setName(text.data());
        
        state_->addVariable(varDihedral);
    }
    else
    {
        printf("Syntax Error: Second argument should be the variable type!");
    }
}
