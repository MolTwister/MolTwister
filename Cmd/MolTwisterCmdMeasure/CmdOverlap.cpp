#include "CmdOverlap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"

std::string CCmdOverlap::getCmd()
{
    return "overlap";
}

std::vector<std::string> CCmdOverlap::getCmdLineKeywords()
{
    return { "overlap", "sel" };
}

std::vector<std::string> CCmdOverlap::getCmdHelpLines()
{
    return {
                "overlap <within> sel"
           };
}

std::string CCmdOverlap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCounts number of overlaping atoms within the current frame by counting the\r\n";
    text+= "\tnumber of atoms that are closer than the distance <within>. The count is\r\n";
    text+= "\tlimited to the currently selected atoms, which is specified by the 'sel'\r\n";
    text+= "\tkeyword.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. Total number of overlaps within <within> = N\r\n";
    text+= "\t2.\r\n";
    text+= "\t3. Atom i  Atom j  Type i  Type j  r_ij\r\n";
    text+= "\t4. ---------------------------------------------------------------------------\r\n";
    text+= "\t5. <index i> <index j> <atom ID i> <atom ID j> <distance from i to j>\r\n";
    text+= "\t6. <index i> <index j> <atom ID i> <atom ID j> <distance from i to j>\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\t       .\r\n";
    text+= "\tN+4. <index i> <index j> <atom ID i> <atom ID j> <distance from i to j>\r\n";
    text+= "\twhere N is the number of identified overlaps.";

    return text;
}

std::string CCmdOverlap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> overlapReport;
    std::string text;
    double withinR2;
    char sstringReport[512];
    bool selectOverlapping = false;
    int overlapCount = 0;

    text = CASCIIUtility::getArg(arguments, arg++);
    withinR2 = atof(text.data());
    withinR2*= withinR2;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "sel")
    {
        selectOverlapping = true;

        for(int i=0; i<state_->atoms_.size(); i++)
        {
            state_->atoms_[i]->select(false);
        }
    }
    else arg--;

    pb.beginProgress("Counting number of atomic overlaps");
    for(int i=0; i<(state_->atoms_.size()-1); i++)
    {
        for(int j=i+1; j<state_->atoms_.size(); j++)
        {
            C3DVector r1 = state_->atoms_[i]->r_[state_->currentFrame_];
            C3DVector r2 = state_->atoms_[j]->r_[state_->currentFrame_];
            C3DVector vDeltaR = r2 - r1;
            double deltaR = vDeltaR.norm2();

            if(deltaR < withinR2)
            {
                overlapCount++;

                std::string ID1 = state_->atoms_[i]->getID();
                std::string ID2 = state_->atoms_[j]->getID();
                sprintf(sstringReport, "\t%-15i%-15i%-15s%-15s%-15g\r\n", i, j, ID1.data(), ID2.data(), sqrt(deltaR));
                overlapReport.emplace_back(sstringReport);

                if(selectOverlapping) state_->atoms_[j]->select();
            }
        }
        pb.updateProgress(i, ((int)state_->atoms_.size()-1));
    }
    pb.endProgress();

    fprintf(stdOut_, "\r\n\tTotal number of overlaps within %g = %i\r\n\r\n", sqrt(withinR2), overlapCount);
    fprintf(stdOut_, "\t%-15s%-15s%-15s%-15s%-15s\r\n\t---------------------------------------------------------------------------\r\n", "Atom i", "Atom j", "Type i", "Type j", "r_ij");
    for(int i=0; i<overlapReport.size(); i++)
    {
        fprintf(stdOut_, "%s", overlapReport[i].data());
    }

    return lastError_;
}
