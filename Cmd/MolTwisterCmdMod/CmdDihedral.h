#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

namespace MolTwisterCmdMod
{
    class CCmdDihedral : public CCmdEntry
    {
    public:
        CCmdDihedral() = delete;
        CCmdDihedral(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
        virtual ~CCmdDihedral() = default;

    public:
        std::string getCmd();
        std::vector<std::string> getCmdLineKeywords();
        std::vector<std::string> getCmdHelpLines();
        std::string getCmdFreetextHelp();
        std::string execute(std::vector<std::string> arguments);

    private:
        std::string lastError_;
    };
}
