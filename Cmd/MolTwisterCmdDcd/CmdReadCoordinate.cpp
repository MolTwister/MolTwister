#include "CmdReadCoordinate.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdReadCoordinate::getCmd()
{
    return "readcoordinate";
}

std::vector<std::string> CCmdReadCoordinate::getCmdLineKeywords()
{
    return { "readcoordinate" };
}

std::vector<std::string> CCmdReadCoordinate::getCmdHelpLines()
{
    return {
                "readcoordinate <DCD filename> <record index> <coordinate index>"
           };
}

std::string CCmdReadCoordinate::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the coordinate at record index <record index>, coordinate index\r\n";
    text+= "\t<coordinate index> within the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tCoordinate = (<x>, <y>, <z>)";

    return text;
}

std::string CCmdReadCoordinate::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    int recordIndex;
    int coordIndex;
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

    text = CASCIIUtility::getArg(arguments, arg++);
    coordIndex = atoi(text.data());

    dcdFile.gotoRecord(recordIndex);

    if(coordIndex >= dcdFile.getNumCoordinatesInRecord())
    {
        lastError_ = std::string("Error coordinate index ") + std::to_string(coordIndex) + std::string(" does not exist!");
        return lastError_;
    }

    C3DVector v = dcdFile.getCoordinate(coordIndex);
    fprintf(stdOut_, "\r\n\tCoordinate = (%.4f, %.4f, %.4f)\r\n", v.x_, v.y_, v.z_);

    return lastError_;
}
