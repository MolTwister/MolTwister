#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

namespace MolTwisterCmdMod
{
    class CCmdAtomPos : public CCmdEntry
    {
    public:
        CCmdAtomPos() = delete;
        CCmdAtomPos(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
        virtual ~CCmdAtomPos() = default;

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
