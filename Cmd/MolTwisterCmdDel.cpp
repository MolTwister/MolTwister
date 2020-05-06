#include <iostream>
#include <vector>
#include "MolTwisterCmdDel.h"

void CCmdDel::onAddKeywords()
{
    addKeyword("del");
    addKeyword("atom");
    addKeyword("atomname");
    addKeyword("sel");
    addKeyword("bond");
    addKeyword("bonds");
    addKeyword("mdnonbonded");
    addKeyword("mdbond");
    addKeyword("mdangle");
    addKeyword("mddihedral");
    addKeyword("ff");
}

std::string CCmdDel::getHelpString() const
{ 
    std::string  text;
    
    text+= "\tUsage: del <specifier>\r\n";
    text+= "\r\n";
    text+= "\tDelete atoms, bonds, angles and dihedrals. The allowed specifiers are:\r\n";
    text+= "\r\n";
    text+= "\t        * atom <N>             :   Delete atom index <N>\r\n";
    text+= "\t        * atomname <name>      :   Delete atoms with name <name>\r\n";
    text+= "\t        * sel                  :   Delete selected atoms\r\n";
    text+= "\t        * bond <N1> <N2>       :   Delete bond <N1>-<N2>\r\n";
    text+= "\t        * bonds <1> <2>        :   Delete bonds <1>-<2>, where <1>, <2> are atom names.\r\n";
    text+= "\t                                   If <2>=* then all bonds to <1> are deleted\r\n";
    text+= "\t        * mdnonbonded <index>  :   The 'list ff' command yields list of indices\r\n";
    text+= "\t        * mdbond <index>       :   The 'list ff' command yields list of indices\r\n";
    text+= "\t        * mdangle <index>      :   The 'list ff' command yields list of indices\r\n";
    text+= "\t        * mddihedral <index>   :   The 'list ff' command yields list of indices\r\n";
    text+= "\t        * ff                   :   Delete force-field (m, Q, mobility, non-bonded, bonds, angles etc.)";
    
    return text;
}

void CCmdDel::execute(std::string commandLine)
{
    std::string text;
    int arg = 1;
    
    if(!state_) return;
    
    text = CASCIIUtility::getWord(commandLine, arg++);
    if(text == "atom")
    {
        int     iAtomIndex;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        iAtomIndex = atoi(text.data());
        if(iAtomIndex < state_->atoms_.size())
        {
            state_->deleteAtom(iAtomIndex);
        }
        else
        {
            printf("Invalid atom index!");
        }
    }
    else if(text == "atomname")
    {
        int size = (int)state_->atoms_.size();
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        for(int i=0; i<size; i++)
        {
            std::string ID = state_->atoms_[i]->getID();
            if(ID == text)
            {
                state_->deleteAtom(i);
                size--; i--;
            }
        }
    }
    else if(text == "sel")
    {
        int size = (int)state_->atoms_.size();
        
        for(int i=0; i<size; i++)
        {
            if(state_->atoms_[i]->isSelected())
            {
                state_->deleteAtom(i);
                size--; i--;
            }
        }
    }
    else if(text == "bond")
    {
        int atomIndex1, atomIndex2;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex1 = atoi(text.data());
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        atomIndex2 = atoi(text.data());
        
        if((atomIndex1 < state_->atoms_.size()) && (atomIndex2 < state_->atoms_.size()))
        {
            CAtom*  pAt1 = state_->atoms_[atomIndex1].get();
            CAtom*  pAt2 = state_->atoms_[atomIndex2].get();
            
            pAt1->deleteBondTo(pAt2);
            pAt2->deleteBondTo(pAt1);
        }
        else
        {
            printf("Invalid atom indices!");
        }
    }
    else if(text == "bonds")
    {
        bool del;
        
        std::string name1 = CASCIIUtility::getWord(commandLine, arg++);
        std::string name2 = CASCIIUtility::getWord(commandLine, arg++);

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            CAtom* at1Ptr = state_->atoms_[i].get();
            std::string atom1 = at1Ptr->getID();
            
            if(atom1 == name1)
            {
                for(int j=0; j<state_->atoms_.size(); j++)
                {
                    del = false;
                    CAtom* at2Ptr = state_->atoms_[j].get();
                    if(at1Ptr == at2Ptr) continue;
                    std::string atom2 = at2Ptr->getID();
                    
                    if(atom2 == name2) del = true;
                    if(name2 == "*") del = true;

                    if(del)
                    {
                        at1Ptr->deleteBondTo(at2Ptr);
                        at2Ptr->deleteBondTo(at1Ptr);
                    }
                }
            }
        }        
    }
    else if(text == "mdnonbonded")
    {
        int index;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        index = atoi(text.data());

        if(index < state_->mdFFNonBondedList_.size())
        {
            state_->mdFFNonBondedList_.del(index);
        }
        else
        {
            printf("Index %i does not exist in list of Non-bonded forcefield parameters!", index);
        }
    }
    else if(text == "mdbond")
    {
        int index;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        index = atoi(text.data());
        
        if(index < state_->mdFFBondList_.size())
        {
            state_->mdFFBondList_.del(index);
        }
        else
        {
            printf("Index %i does not exist in list of bond forcefield parameters!", index);
        }
    }
    else if(text == "mdangle")
    {
        int index;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        index = atoi(text.data());
        
        if(index < state_->mdFFAngleList_.size())
        {
            state_->mdFFAngleList_.del(index);
        }
        else
        {
            printf("Index %i does not exist in list of angle forcefield parameters!", index);
        }
    }
    else if(text == "mddihedral")
    {
        int index;
        
        text = CASCIIUtility::getWord(commandLine, arg++);
        index = atoi(text.data());
        
        if(index < state_->mdFFDihList_.size())
        {
            state_->mdFFDihList_.del(index);
        }
        else
        {
            printf("Index %i does not exist in list of dihedral forcefield parameters!", index);
        }
    }
    else if(text == "ff")
    {
        state_->mdFFNonBondedList_.empty();
        state_->mdFFBondList_.empty();
        state_->mdFFAngleList_.empty();
        state_->mdFFDihList_.empty();
        
        for(int i=0; i<state_->atoms_.size(); i++)
        {
            state_->atoms_[i]->isMobile_ = true;
            state_->atoms_[i]->Q_ = 0.0;
            state_->atoms_[i]->m_ = 0.0;
        }
    }
    else
    {
        printf("Syntax Error: First argument should be the type of object to delete!");
    }
    
    if(state_->view3D_)
    {
        state_->view3D_->requestUpdate(false);
    }
    
    for(int i=0; i<state_->atoms_.size(); i++)
    {
        if(state_->atoms_[i])
            state_->atoms_[i]->buildListOf1to4Connections();
    }
}
