#include "CmdWrap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../../Utilities/FileUtility.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdWrap::getCmd()
{
    return "wrap";
}

std::vector<std::string> CCmdWrap::getCmdLineKeywords()
{
    return { "wrap", "pbcfromvisual" };
}

std::vector<std::string> CCmdWrap::getCmdHelpLines()
{
    return {
                "wrap <DCD filename> [pbcfromvisual]"
           };
}

std::string CCmdWrap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tLoads a DCD file, <DCD filename>, and performs a wrap across the periodic boundaries\r\n";
    text+= "\t(i.e., PBC). The PBC is taken from the DCD file by default. However, by specifying\r\n";
    text+= "\tthe 'pbcfromvisual' keyword, the PBC can be taken from the applied visual PBC.\r\n";
    text+= "\tThe output is a DCD file with the name\r\n";
    text+= "\t* <DCD filename (without extension)>_mtwrap.dcd";

    return text;
}

std::string CCmdWrap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::string text;
    C3DRect pbc;
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
    outFileName+= "_mtwrap.dcd";

    FILE* outFile = fopen(outFileName.data(), "w");
    if(!outFile)
    {
        lastError_ = std::string("Error opening file ") + outFileName + std::string(" for writing!");
        return lastError_;
    }


    // Write a copy of the main header to output dcd. If requested
    // choose to use PBC displayed in the MolTwister 3D view, instead of PBC
    // from dcd file record headers (i.e. for use in wrapping algorithm)
    int numBytesWritten;
    bool usePBCFromVisual = false;
    CDCDFile::CMainHeader Header = dcdFile.getDCDHeader();
    Header.write(outFile, numBytesWritten);

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "pbcfromvisual")
    {
        pbc = state_->view3D_->calcPBC();
        usePBCFromVisual = true;

    } else arg--;


    // Check that we have any records and not the maximum
    // number of coordinates (taken from first record)
    int numRecords = dcdFile.getNumRecords();
    if(numRecords <= 0)
    {
        fclose(outFile);
        lastError_ = "Error: found no records!";
        return lastError_;
    }

    dcdFile.gotoRecord(0);
    int maxCoordinates = dcdFile.getNumCoordinatesInRecord();


    // Wrap coordinates and store to output dcd file
    C3DVector r;
    pb.beginProgress("Wrapping DCD file");
    for(int t=0; t<numRecords; t++)
    {
        dcdFile.gotoRecord(t);

        if(!usePBCFromVisual)
        {
            pbc.rHigh_.x_ = dcdFile.getCurrentRecordHeader().boxX_ / 2.0;
            pbc.rHigh_.y_ = dcdFile.getCurrentRecordHeader().boxY_ / 2.0;
            pbc.rHigh_.z_ = dcdFile.getCurrentRecordHeader().boxZ_ / 2.0;

            pbc.rLow_.x_ = -pbc.rHigh_.x_;
            pbc.rLow_.y_ = -pbc.rHigh_.y_;
            pbc.rLow_.z_ = -pbc.rHigh_.z_;
        }

        int iNumCoordinates = dcdFile.getNumCoordinatesInRecord();
        for(int i=0; i<iNumCoordinates; i++)
        {
            // Read current coordinate and calculate
            // wrapped coordinates
            if(i >= maxCoordinates) continue;
            r = dcdFile.getCoordinate(i);

            C3DVector add;
            if(r.x_ > pbc.rHigh_.x_) add.x_ = -double(int((r.x_ - pbc.rLow_.x_) / pbc.getWidthX())) * pbc.getWidthX();
            if(r.x_ < pbc.rLow_.x_)  add.x_ =  double(int((pbc.rHigh_.x_ - r.x_) / pbc.getWidthX())) * pbc.getWidthX();
            if(r.y_ > pbc.rHigh_.y_) add.y_ = -double(int((r.y_ - pbc.rLow_.y_) / pbc.getWidthY())) * pbc.getWidthY();
            if(r.y_ < pbc.rLow_.y_)  add.y_ =  double(int((pbc.rHigh_.y_ - r.y_) / pbc.getWidthY())) * pbc.getWidthY();
            if(r.z_ > pbc.rHigh_.z_) add.z_ = -double(int((r.z_ - pbc.rLow_.z_) / pbc.getWidthZ())) * pbc.getWidthZ();
            if(r.z_ < pbc.rLow_.z_)  add.z_ =  double(int((pbc.rHigh_.z_ - r.z_) / pbc.getWidthZ())) * pbc.getWidthZ();

            // Update current coordinate to wrapped state
            dcdFile.setCoordinate(i, r + add);
        }

        dcdFile.saveCurrentRecord(outFile);
        pb.updateProgress(t, numRecords);
    }
    pb.endProgress();

    fclose(outFile);

    return lastError_;
}
