#include <iostream>
#include <vector>
#include "MolTwisterCmdVarlist.h"

void CCmdVarlist::onAddKeywords()
{
    addKeyword("varlist");
}

std::string CCmdVarlist::getHelpString() const
{
    std::string text;
    
    text+= "\tUsage: varlist\r\n";
    text+= "\r\n";
    text+= "\tList all variables that have been defined using the 'var command.\r\n";
    text+= "\tThe variable names, their type and their content will be listed.";
    
    return text;
}

void CCmdVarlist::execute(std::string)
{
    std::string text;
    
    if(!state_) return;
    
    fprintf(stdOut_, "\r\n\r\n\t%-40s%-20s%-10s\r\n", "Variable", "Type", "Contents");
    fprintf(stdOut_, "\t----------------------------------------------------------------------------------------\r\n");
    for(int i=0; i<state_->variables_.size(); i++)
    {
        if(state_->variables_[i]->getType() == CVar::typeAtom)
        {
            CVarAtom* p = (CVarAtom*)state_->variables_[i].get();
            
            state_->variables_[i]->getName(text);
            fprintf(stdOut_, "\t%-40s%-20sAtom Index=%i\r\n", text.data(), "[atom]", p->atomIndex_);
        }
        if(state_->variables_[i]->getType() == CVar::typeBond)
        {
            CVarBond* p = (CVarBond*)state_->variables_[i].get();
            
            state_->variables_[i]->getName(text);
            fprintf(stdOut_, "\t%-40s%-20sBond Atom Indices={%i,%i}\r\n", text.data(), "[bond]", p->atomIndex1_, p->atomIndex2_);
        }
        if(state_->variables_[i]->getType() == CVar::typeAngle)
        {
            CVarAngle* p = (CVarAngle*)state_->variables_[i].get();
            
            state_->variables_[i]->getName(text);
            fprintf(stdOut_, "\t%-40s%-20sAngle Atom Indices={%i,%i,%i}\r\n", text.data(), "[angle]", p->atomIndex1_, p->atomIndex2_, p->atomIndex3_);
        }
        if(state_->variables_[i]->getType() == CVar::typeDihedral)
        {
            CVarDihedral* p = (CVarDihedral*)state_->variables_[i].get();
            
            state_->variables_[i]->getName(text);
            fprintf(stdOut_, "\t%-40s%-20sDihedral Atom Indices={%i,%i,%i,%i}\r\n", text.data(), "[dihedral]", p->atomIndex1_, p->atomIndex2_, p->atomIndex3_, p->atomIndex4_);
        }
    }
    fprintf(stdOut_, "\r\n");
}
