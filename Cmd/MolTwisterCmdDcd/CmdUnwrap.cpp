#include "CmdUnwrap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../../Utilities/FileUtility.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdUnwrap::getCmd()
{
    return "unwrap";
}

std::vector<std::string> CCmdUnwrap::getCmdLineKeywords()
{
    return { "unwrap", "pbcfromvisual" };
}

std::vector<std::string> CCmdUnwrap::getCmdHelpLines()
{
    return {
                "unwrap <DCD filename> [pbcfromvisual]"
           };
}

std::string CCmdUnwrap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tLoads a DCD file, <DCD filename>, and performs an unwrap across the periodic\r\n";
    text+= "\tboundaries (i.e., PBC). The PBC is taken from the DCD file by default. However,\r\n";
    text+= "\tby specifying the 'pbcfromvisual' keyword, the PBC can be taken from the applied\r\n";
    text+= "\tvisual PBC. The output is a DCD file with the name\r\n";
    text+= "\t* <DCD filename (without extension)>_mtunwrap.dcd\r\n";
    text+= "\r\n";
    text+= "\tNote that an 'unwrap' will move the atomic coordinates across the PBC if it\r\n";
    text+= "\tis discovered that the coordinate crosses it between two time steps (thus leading\r\n";
    text+= "\tto a complete unfolding of the PBC. However, an 'atomic unwrap' will unwrap\r\n";
    text+= "\tframe-by-frame, based on the molecular definitions. Hence, only molecules that\r\n";
    text+= "\twrap across the PBC will be unwrapped so as to be un-divided, but outside the\r\n";
    text+= "\tPBC. Hence. an atomic unwrap requires only one frame, while unwrap needs more.";

    return text;
}

std::string CCmdUnwrap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    double pbc[3];
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
    outFileName+= "_mtunwrap.dcd";

    FILE* outFile = fopen(outFileName.data(), "w");
    if(!outFile)
    {
        lastError_ = std::string("Error opening file ") + outFileName + std::string(" for writing!");
        return lastError_;
    }


    // Write a copy of the main header to output dcd. If requested
    // choose to use PBC displayed in the MolTwister 3D view, instead of PBC
    // from dcd file record headers (i.e. for use in unwrapping algorithm)
    int numBytesWritten;
    bool usePBCFromVisual = false;
    CDCDFile::CMainHeader header = dcdFile.getDCDHeader();
    header.write(outFile, numBytesWritten);

    C3DRect visualPBC;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "pbcfromvisual")
    {
        visualPBC = state_->view3D_->calcPBC();
        pbc[0] = visualPBC.getWidthX();
        pbc[1] = visualPBC.getWidthY();
        pbc[2] = visualPBC.getWidthZ();
        usePBCFromVisual = true;

    } else arg--;


    // Prepare vector keeping track of the PBC jumps of all
    // the coordinates in a system, as well as the array
    // containing last coordinates
    int numRecords = dcdFile.getNumRecords();
    if(numRecords <= 0)
    {
        fclose(outFile);
        lastError_ = "Error: found no records!";
        return lastError_;
    }

    dcdFile.gotoRecord(0);
    int maxCoordinates = dcdFile.getNumCoordinatesInRecord();
    std::vector<C3DVector> pbcJumpDistances;
    std::vector<C3DVector> lastCoordinates;
    pbcJumpDistances.resize(maxCoordinates);
    lastCoordinates.resize(maxCoordinates);
    for(int i=0; i<maxCoordinates; i++) { lastCoordinates[i] = dcdFile.getCoordinate(i); }


    // Unwrap coordinates and store to output dcd file
    C3DVector r1, r2;
    pb.beginProgress("Unwrapping DCD file");
    for(int t=0; t<numRecords; t++)
    {
        dcdFile.gotoRecord(t);

        if(!usePBCFromVisual)
        {
            pbc[0] = dcdFile.getCurrentRecordHeader().boxX_;
            pbc[1] = dcdFile.getCurrentRecordHeader().boxY_;
            pbc[2] = dcdFile.getCurrentRecordHeader().boxZ_;
        }

        int numCoordinates = dcdFile.getNumCoordinatesInRecord();
        for(int i=0; i<numCoordinates; i++)
        {
            // Read current and previous coordinate,
            // check if it changes PBC, if so
            // update jump distance accordingly.
            // Also, store last coordinate before update
            if(i >= maxCoordinates) continue;
            r1 = lastCoordinates[i];
            r2 = dcdFile.getCoordinate(i);

            double dx = r2.x_ - r1.x_;
            double dy = r2.y_ - r1.y_;
            double dz = r2.z_ - r1.z_;

            C3DVector add;
            if(fabs(dx) > (pbc[0]/2.0))
            {
                if(dx < 0.0) add.x_ = pbc[0];
                else         add.x_ = -pbc[0];
            }
            if(fabs(dy) > (pbc[1]/2.0))
            {
                if(dy < 0.0) add.y_ = pbc[1];
                else         add.y_ = -pbc[1];
            }
            if(fabs(dz) > (pbc[2]/2.0))
            {
                if(dz < 0.0) add.z_ = pbc[2];
                else         add.z_ = -pbc[2];
            }

            pbcJumpDistances[i].x_+= add.x_;
            pbcJumpDistances[i].y_+= add.y_;
            pbcJumpDistances[i].z_+= add.z_;

            lastCoordinates[i] = r2;

            // Update current coordinate to unwrapped state
            dcdFile.setCoordinate(i, r2 + pbcJumpDistances[i]);
        }

        dcdFile.saveCurrentRecord(outFile);
        pb.updateProgress(t, numRecords);
    }
    pb.endProgress();

    fclose(outFile);

    return lastError_;
}
