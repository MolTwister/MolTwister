#include "CmdNumRecords.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdNumRecords::getCmd()
{
    return "numrecords";
}

std::vector<std::string> CCmdNumRecords::getCmdLineKeywords()
{
    return { "numrecords" };
}

std::vector<std::string> CCmdNumRecords::getCmdHelpLines()
{
    return {
                "numrecords <DCD filename>"
           };
}

std::string CCmdNumRecords::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the number of recoords within the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tNum records = <record count>";

    return text;
}

std::string CCmdNumRecords::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CDCDFile dcdFile;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error opening file ") + text + std::string("!");
        return lastError_;
    }

    fprintf(stdOut_, "\r\n\tNum records = %i\r\n", dcdFile.getNumRecords());

    return lastError_;
}
