#include "CmdNumCoordinates.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdNumCoordinates::getCmd()
{
    return "numcoordinates";
}

std::vector<std::string> CCmdNumCoordinates::getCmdLineKeywords()
{
    return { "numcoordinates" };
}

std::vector<std::string> CCmdNumCoordinates::getCmdHelpLines()
{
    return {
                "numcoordinates <DCD filename> <record index>"
           };
}

std::string CCmdNumCoordinates::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the number of coordinates in record <record index>\r\n";
    text+= "\tof the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tCoordinate count = <coordinate count>";

    return text;
}

std::string CCmdNumCoordinates::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    int recordIndex;
    std::string text;
    CDCDFile dcdFile;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error opening file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    recordIndex = atoi(text.data());

    dcdFile.gotoRecord(recordIndex);

    fprintf(stdOut_, "\r\n\tCoordinate count = %i\r\n", dcdFile.getNumCoordinatesInRecord());

    return lastError_;
}
