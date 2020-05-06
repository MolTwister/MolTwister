#include "CmdAtomicUnwrap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../../Utilities/FileUtility.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdAtomicUnwrap::getCmd()
{
    return "atomicunwrap";
}

std::vector<std::string> CCmdAtomicUnwrap::getCmdLineKeywords()
{
    return { "atomicunwrap", "atomnames", "pbcfromvisual" };
}

std::vector<std::string> CCmdAtomicUnwrap::getCmdHelpLines()
{
    return {
                "atomicunwrap <DCD filename> [atomnames <atom ID list>] [pbcfromvisual]"
           };
}

std::string CCmdAtomicUnwrap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tLoads a DCD file, <DCD filename>, and performs an atomic unwrap across the\r\n";
    text+= "\tperiodic boundaries (i.e., PBC). The PBC is taken from the DCD file by defualt.\r\n";
    text+= "\tHowever, by specifying the 'pbcfromvisual' keyword, the PBC can be taken from\r\n";
    text+= "\tthe applied visual PBC. It is also possible to select a subsection of the atoms\r\n";
    text+= "\tstored in the original DCD file by applying the 'atomnames' keyword with a list\r\n";
    text+= "\tof atom IDs (e.g., H, O, C7), <atom ID list>, which is comma separated with no\r\n";
    text+= "\tspace. The output is a DCD file with the name\r\n";
    text+= "\t* <DCD filename (without extension)>_mtatunwrap.dcd\r\n";
    text+= "\r\n";
    text+= "\tNote that an 'unwrap' will move the atomic coordinates across the PBC if it\r\n";
    text+= "\tis discovered that the coordinate crosses it between two time steps (thus leading\r\n";
    text+= "\tto a complete unfolding of the PBC. However, an 'atomic unwrap' will unwrap\r\n";
    text+= "\tframe-by-frame, based on the molecular definitions. Hence, only molecules that\r\n";
    text+= "\twrap across the PBC will be unwrapped so as to be un-divided, but outside the\r\n";
    text+= "\tPBC. Hence. an atomic unwrap requires only one frame, while unwrap needs more.";

    return text;
}

std::string CCmdAtomicUnwrap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::string text;
    CDCDFile dcdFile;


    // Open original dcd file
    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error opening file ") + text + std::string("!");
        return lastError_;
    }


    // Construct output dcd file name and create the file
    std::string outFileName = text;
    CFileUtility::removeExtension(outFileName);
    outFileName+= "_mtatunwrap.dcd";

    FILE* outFile = fopen(outFileName.data(), "w");
    if(!outFile)
    {
        lastError_ = std::string("Error opening file ") + outFileName + std::string(" for writing!");
        return lastError_;
    }


    // Retrieve a list of atom types to consider (or all, if not specifically chosen)
    // and use this list to retrieve an appropriate molecular list
    std::vector<std::string> atomTypesToInclude;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "atomnames")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        CASCIIUtility::removeWhiteSpace(text);
        atomTypesToInclude = CASCIIUtility::getWords(text, ",");

    } else arg--;

    std::vector<std::vector<int>> molList;
    CMolTwisterStateTools(state_, stdOut_).getMolecules(molList, atomTypesToInclude);



    // Write a copy of the main header to output dcd. If requested
    // choose to use PBC displayed in the MolTwister 3D view, instead of PBC
    // from dcd file record headers (i.e. for use in unwrapping algorithm)
    int numBytesWritten;
    bool usePBCFromVisual = false;
    CDCDFile::CMainHeader Header = dcdFile.getDCDHeader();
    Header.write(outFile, numBytesWritten);

    C3DRect visualPBC;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "pbcfromvisual")
    {
        visualPBC = state_->view3D_->calcPBC();
        usePBCFromVisual = true;

    } else arg--;


    // Get number of records in DCD file
    int numRecords = dcdFile.getNumRecords();
    if(numRecords <= 0)
    {
        fclose(outFile);
        lastError_ = "Error: found no records!";
        return lastError_;
    }


    // Unwrap coordinates and store to output dcd file
    CDCDTools dcdTools(state_, stdOut_);
    pb.beginProgress("Unwrapping DCD file");
    for(int t=0; t<numRecords; t++)
    {
        dcdFile.gotoRecord(t);

        if(usePBCFromVisual) dcdTools.atomicUnwrapCurrDCDFrame(&dcdFile, molList, &visualPBC);
        else                 dcdTools.atomicUnwrapCurrDCDFrame(&dcdFile, molList, nullptr);

        dcdFile.saveCurrentRecord(outFile);
        pb.updateProgress(t, numRecords);
    }
    pb.endProgress();

    fclose(outFile);

    return lastError_;
}
