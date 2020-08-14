#pragma once
#include "../Tools/MolTwisterCmdEntry.h"

namespace MolTwisterCmdMeasure
{
    class CCmdAngle : public CCmdEntry
    {
    public:
        CCmdAngle() = delete;
        CCmdAngle(CMolTwisterState* state, FILE* stdOut) : CCmdEntry(state, stdOut) { }
        virtual ~CCmdAngle() = default;

    public:
        std::string getCmd();
        std::vector<std::string> getCmdLineKeywords();
        std::vector<std::string> getCmdHelpLines();
        std::string getCmdFreetextHelp();
        std::string execute(std::vector<std::string> arguments);
        void redirectOutput(FILE* stdOut = stdout) { stdOut_ = stdOut; }

    private:
        std::string lastError_;
    };
}
