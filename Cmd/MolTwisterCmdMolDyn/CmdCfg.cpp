#include "CmdCfg.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdCfg::getCmd()
{
    return "cfg";
}

std::vector<std::string> CCmdCfg::getCmdLineKeywords()
{
    std::vector<std::string> keywords = { "cfg", "get", "set", "resettodefaults" };
    std::vector<std::string> keywordsCfg = molDynConfig_->getKeyWords();

    for(auto keyword : keywordsCfg)
    {
        keywords.emplace_back(keyword);
    }

    return keywords;
}

std::vector<std::string> CCmdCfg::getCmdHelpLines()
{
    return {
        "cfg get",
        "cfg set <parameter> <value>"
    };
}

std::string CCmdCfg::getCmdFreetextHelp()
{
    std::string text;

    text+= "\t\r\n";

    return text;
}

std::string CCmdCfg::execute(std::vector<std::string> arguments)
{
    lastError_ = "";
    size_t arg = 0;

    std::string cmd =  CASCIIUtility::getArg(arguments, arg++).data();
    if(cmd == "get")
    {
        parseGetCommand(arguments, arg);
    }
    else if(cmd == "set")
    {
        parseSetCommand(arguments, arg);
    }
    else if(cmd == "resettodefaults")
    {
        parseResettodefaultsCommand(arguments, arg);
    }

    return lastError_;
}

void CCmdCfg::parseGetCommand(const std::vector<std::string>&, size_t&)
{
    molDynConfig_->print(stdOut_);
}

void CCmdCfg::parseSetCommand(const std::vector<std::string>& arguments, size_t& arg)
{
    std::string parameter = CASCIIUtility::getArg(arguments, arg++);
    std::string value = CASCIIUtility::getArg(arguments, arg++);

    lastError_ = molDynConfig_->set(parameter, value);
}

void CCmdCfg::parseResettodefaultsCommand(const std::vector<std::string>&, size_t&)
{
    molDynConfig_->resetToDefaults();
}