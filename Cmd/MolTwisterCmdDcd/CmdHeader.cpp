#include "CmdHeader.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"

std::string CCmdHeader::getCmd()
{
    return "header";
}

std::vector<std::string> CCmdHeader::getCmdLineKeywords()
{
    return { "header" };
}

std::vector<std::string> CCmdHeader::getCmdHelpLines()
{
    return {
                "header <DCD filename>"
           };
}

std::string CCmdHeader::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tPrints the DCD file header from <DCD filename>.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\tID = <ID>\r\n";
    text+= "\tNSets = <number of sets>\r\n";
    text+= "\tInitStep = <init step>\r\n";
    text+= "\tWrtFreq = <Write frequency>\r\n";
    text+= "\tTimeStep = <time step>\r\n";
    text+= "\tDescriptA = <descript A>\r\n";
    text+= "\tDescriptB = <descript B>\r\n";
    text+= "\tNAtoms = <number of atoms>";

    return text;
}

std::string CCmdHeader::execute(std::vector<std::string> arguments)
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

    fprintf(stdOut_, "\r\n");
    fprintf(stdOut_, "\tID = %s\r\n", dcdFile.getDCDHeader().ID_.data());
    fprintf(stdOut_, "\tNSets = %i\r\n", dcdFile.getDCDHeader().nSets_);
    fprintf(stdOut_, "\tInitStep = %i\r\n", dcdFile.getDCDHeader().initStep_);
    fprintf(stdOut_, "\tWrtFreq = %i\r\n", dcdFile.getDCDHeader().wrtFreq_);
    fprintf(stdOut_, "\tTimeStep = %g\r\n", (double)dcdFile.getDCDHeader().timeStep_);
    fprintf(stdOut_, "\tDescriptA = %s\r\n", dcdFile.getDCDHeader().descriptA_.data());
    fprintf(stdOut_, "\tDescriptB = %s\r\n", dcdFile.getDCDHeader().descriptB_.data());
    fprintf(stdOut_, "\tNAtoms = %i\r\n", dcdFile.getDCDHeader().nAtoms_);

    return lastError_;
}
