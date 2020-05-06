#include "CmdMDDihedral.h"
#include "../../Utilities/ASCIIUtility.h"

#define DEFAULT_MAX_BOND 1.75

std::string CCmdMDDihedral::getCmd()
{
    return "mddihedral";
}

std::vector<std::string> CCmdMDDihedral::getCmdLineKeywords()
{
    return { "mddihedral" };
}

std::vector<std::string> CCmdMDDihedral::getCmdHelpLines()
{
    return {
                "mddihedral <ID1> <ID2> <ID3> <ID4> <FF-type> [all_r_less_than <radius>, mol_r_less_than <radius>, only_visible_bonds] <parameters for given FF-type>"
           };
}

std::string CCmdMDDihedral::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tID1, ID2, ID3 and ID4 identifiy the atom types where a dihedral interaction is to be\r\n";
    text+= "\tdefined (e.g., H, O, C, O5). The force field type, FF-type, is idefined by\r\n";
    text+= "\ta string. The possible strings are listed below, together with how 'parameters\r\n";
    text+= "\tfor given FF-type' is defined for eaah string.\r\n";
    text+= "\r\n";
    text+= "\tPossible bond definitions to apply:\r\n";
    text+= "\t* all_r_less_than <r> - yields bonds even between atoms not close enough to define a molecule (if the r-criteria is satisfied)\r\n";
    text+= "\t* mol_r_less_than <r> - yields bonds only between atoms close enough to define a molecule (if the r-criteria is satisfied)\r\n";
    text+= "\t* only_visible_bonds - yields bonds only where bonds are visible in the 3D view\r\n";

    text+= "\r\n";

    text+= "\tPossible <FF-type> <parameters for given FF-type> to apply:\r\n";

    int regTypes = state_->mdFFDihList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        CMDFFDih* ffType = state_->mdFFDihList_.getRegisteredFFType(i);
        if(ffType)
        {
            std::string ffTypeString = ffType->getFFType();
            std::vector<std::string> paramsHelpLines = ffType->getCmdHelpLines();

            for(std::string paramsHelp : paramsHelpLines)
            {
                text+= std::string("\t* ") + ffTypeString + std::string(" ") + paramsHelp + std::string("\r\n");
            }
        }
    }

    return text;
}

std::string CCmdMDDihedral::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CMDFFDih::EDetCrit detCrit;
    std::shared_ptr<CMDFFDih> dihEntry;
    std::string stringAtom1, stringAtom2, stringAtom3, stringAtom4, type, text;
    double detCritR0;
    int regTypes;

    stringAtom1 = CASCIIUtility::getArg(arguments, arg++);
    stringAtom2 = CASCIIUtility::getArg(arguments, arg++);
    stringAtom3 = CASCIIUtility::getArg(arguments, arg++);
    stringAtom4 = CASCIIUtility::getArg(arguments, arg++);
    type = CASCIIUtility::getArg(arguments, arg++);

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "all_r_less_than")
    {
        detCrit = CMDFFDih::critAllRLessThan_R0;
        text = CASCIIUtility::getArg(arguments, arg++);
        detCritR0 = atof(text.data());
    }
    else if(text == "mol_r_less_than")
    {
        detCrit = CMDFFDih::critMolRLessThan_R0;
        text = CASCIIUtility::getArg(arguments, arg++);
        detCritR0 = atof(text.data());
    }
    else if(text == "only_visible_bonds")
    {
        detCrit = CMDFFDih::critOnlyVisibleBonds;
        detCritR0 = 0.0;
    }
    else
    {
        detCrit = CMDFFDih::critAllRLessThan_R0;
        detCritR0 = DEFAULT_MAX_BOND;
        arg--;
    }

    regTypes = state_->mdFFDihList_.getNumRegisteredFFTypes();
    for(int i=0; i<regTypes; i++)
    {
        if(type == state_->mdFFDihList_.getRegisteredFFType(i)->getFFType())
        {
            dihEntry = state_->mdFFDihList_.getRegisteredFFType(i)->createCopy();
            break;
        }
    }

    if(dihEntry)
    {
        dihEntry->setAtomsToBond(stringAtom1, stringAtom2, stringAtom3, stringAtom4);
        dihEntry->setBondDetectionCriteria(detCrit, detCritR0);

        arguments.erase(arguments.begin(), arguments.begin() + arg);
        dihEntry->parse(arguments);
        state_->mdFFDihList_.add(*dihEntry);
    }
    else
    {
        lastError_ = std::string("Error: Did not recognize dihedral type ") + type + std::string("!");
        return lastError_;
    }

    return lastError_;
}
