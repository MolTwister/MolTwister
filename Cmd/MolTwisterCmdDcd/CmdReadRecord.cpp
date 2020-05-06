#include "CmdReadRecord.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdReadRecord::getCmd()
{
    return "readrecord";
}

std::vector<std::string> CCmdReadRecord::getCmdLineKeywords()
{
    return { "readrecord" };
}

std::vector<std::string> CCmdReadRecord::getCmdHelpLines()
{
    return {
                "readrecord <DCD filename> <record index>"
           };
}

std::string CCmdReadRecord::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints all the coordinates at record index <record index>\r\n";
    text+= "\twithin the DCD file, <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t{{<x1>, <y1>, <z1>}, {<x2>, <y2>, <z2>}, ..., {<xN>, <yN>, <zN>}}\r\n";
    text+= "\twhere N is the number of coordinates within the given record.";

    return text;
}

std::string CCmdReadRecord::execute(std::vector<std::string> arguments)
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

    fprintf(stdOut_, "\r\n{");
    for(int i=0; i<dcdFile.getNumCoordinatesInRecord(); i++)
    {
        C3DVector v;

        v = dcdFile.getCoordinate(i);
        if(i != 0) fprintf(stdOut_, ",");
        fprintf(stdOut_, "{%.4f,%.4f,%.4f}", v.x_, v.y_, v.z_);
    }
    fprintf(stdOut_, "}");

    return lastError_;
}
