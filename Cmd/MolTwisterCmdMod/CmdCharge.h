#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

namespace MolTwisterCmdMod
{
    class CCmdCharge : public CCmdEntry
    {
    public:
        CCmdCharge() = delete;
        CCmdCharge(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
        virtual ~CCmdCharge() = default;

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
