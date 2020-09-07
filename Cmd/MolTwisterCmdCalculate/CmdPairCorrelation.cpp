#include "CmdPairCorrelation.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"
#include <float.h>

std::string CCmdPairCorrelation::getCmd()
{
    return "paircorrelation";
}

std::vector<std::string> CCmdPairCorrelation::getCmdLineKeywords()
{
    return { "paircorrelation" };
}

std::vector<std::string> CCmdPairCorrelation::getCmdHelpLines()
{
    return {
                "paircorrelation <DCD filename> <frame from> <frame to> <atom ID 1> <atom ID 2> <num bins> <min. dist> <max. dist>"
           };
}

std::string CCmdPairCorrelation::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the pair correlation function, averaged over frames <frame from> to <frame to>,\r\n";
    text+= "\ttaken from the DCD file <DCD filename>. Pair correlation is calculated between atoms with\r\n";
    text+= "\tID <atom ID 1> and <atom ID 2> (e.g., H, O, C7). The pair correlation function is calculated\r\n";
    text+= "\tfor distances between <min. dist> and <max. dist>, divided into <num bins>.\r\n";
    text+= "\r\n";
    text+= "\tOuptut:\r\n";
    text+= "\t1. From frame: <frame from>\r\n";
    text+= "\t2. To frame: <frame to>\r\n";
    text+= "\t3. Range: [<min. dist>, <max. dist>]\r\n";
    text+= "\t4. Num. bins: <num bins>\r\n";
    text+= "\t5. Distance Count PairCorr RDF IdGas.RDF\r\n";
    text+= "\t6. <distance> <num. atoms in bin> <pair correlation value> <radial distribution function (RDF) value> <ideal gas RDF value>\r\n";
    text+= "\t7. <distance> <num. atoms in bin> <pair correlation value> <radial distribution function (RDF) value> <ideal gas RDF value>\r\n";
    text+= "\t        .\r\n";
    text+= "\t        .\r\n";
    text+= "\t        .\r\n";
    text+= "\tN+5. <distance> <num. atoms in bin> <pair correlation value> <radial distribution function (RDF) value> <ideal gas RDF value>\r\n";
    text+= "\twhere N is the number of applied bins.";

    return text;
}

std::string CCmdPairCorrelation::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<int> profile;
    std::vector<double> pairCorr;
    std::vector<double> rdf;
    std::vector<double> dist;
    std::vector<double> idGasPairCorr;
    std::vector<int> atomIndices1;
    std::vector<int> atomIndices2;
    CDCDFile dcdFile;
    std::string text;
    int frameFrom, frameTo;


    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);
    if(!dcdFile.open(text))
    {
        lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
        return lastError_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    frameFrom = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    frameTo = atoi(text.data());

    if(frameTo <= frameFrom)
    {
        lastError_ = "Error: frame to must be larger than frame from!";
        return lastError_;
    }


    // Obtain list of all atom indices to include in bond length spectrum based on atom IDs
    std::string ID1 = CASCIIUtility::getArg(arguments, arg++);
    std::string ID2 = CASCIIUtility::getArg(arguments, arg++);

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        std::string cmp = state_->atoms_[i]->getID();
        if(cmp == ID1) atomIndices1.emplace_back(i);
        if(cmp == ID2) atomIndices2.emplace_back(i);
    }


    // Obtain number of profile bins and profile boundaries
    text = CASCIIUtility::getArg(arguments, arg++);
    int size = atoi(text.data());
    if(size <= 0)
    {
        lastError_ = "Error: cannot have zero or fewer bins in profile!";
        return lastError_;
    }
    profile.resize(size, 0);
    pairCorr.resize(size, 0.0);
    rdf.resize(size, 0.0);
    dist.resize(size, 0.0);
    idGasPairCorr.resize(size, 0.0);

    text = CASCIIUtility::getArg(arguments, arg++);
    double minDist = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double maxDist = atof(text.data());

    if(fabs(double(minDist - maxDist)) < double(FLT_MIN))
    {
        lastError_ = "Error: cannot have equal minimum and maximum distances in profile!";
        return lastError_;
    }


    // Go through all frames and build profile
    double avgNumDens = 0.0;
    double qvgVol = 0.0;
    int count = 0;
    int numFrames = dcdFile.getNumRecords();
    pb.beginProgress("Calculate pair correlation function");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        C3DRect pbc = dcdFile.getCurrentPBC();
        double V = pbc.getVolume();
        double N = double(atomIndices1.size() + atomIndices2.size());
        if(V == 0.0)
        {
            lastError_ = std::string("Error: encountered zero volume at frame ") + std::to_string(t) + std::string(" of ") + std::to_string(numFrames) + std::string("!");
            return lastError_;
        }
        avgNumDens+= N / V;
        qvgVol+= V;
        count++;

        for(int i=0; i<atomIndices1.size(); i++)
        {
            for(int j=0; j<atomIndices2.size(); j++)
            {
                C3DVector r1 = dcdFile.getCoordinate(atomIndices1[i]);
                C3DVector r2 = dcdFile.getCoordinate(atomIndices2[j]);
                C3DVector r12 = r2 - r1;
                double R12 = r12.norm();

                int n = int((R12 - minDist) / (maxDist - minDist) * double(size));
                if((n >= 0) && (n < size))
                    profile[n]++;
            }
        }

        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();

    if(count == 0)
    {
        lastError_ = "Error: found no valid frames!";
        return lastError_;
    }
    avgNumDens/= double(count);
    qvgVol/= double(count);


    // Create pair correlation and radial distrubution functions
    for(int i=0; i<profile.size(); i++)
    {
        dist[i] = 0.5 * double(2*i+1) * (maxDist - minDist) / double(size) + minDist;
        idGasPairCorr[i] = 4.0 * M_PI * dist[i]*dist[i] * avgNumDens;

        pairCorr[i] = double(profile[i]) / qvgVol; // [#/AA^3]
        rdf[i] = pairCorr[i] / idGasPairCorr[i]; // Unitless
    }


    // Print results
    printf("\r\n");
    fprintf(stdOut_, "From frame: %i\r\n", frameFrom);
    fprintf(stdOut_, "To frame: %i\r\n", frameTo);
    fprintf(stdOut_, "Pair: %s-%s\r\n", ID1.data(), ID2.data());
    fprintf(stdOut_, "Range: [%.4f, %.4f]\r\n", minDist, maxDist);
    fprintf(stdOut_, "Num. bins: %i\r\n\r\n", size);
    fprintf(stdOut_, "%-15s%-15s%-15s%-15s%-15s\r\n", "Distance", "Count", "PairCorr", "RDF", "IdGas.RDF");
    fprintf(stdOut_, "------------------------------\r\n");

    for(int i=0; i<profile.size(); i++)
    {
        fprintf(stdOut_, "%-15.4f%-15i%-15.6f%-15.6f%-15.6f\r\n", dist[i], profile[i], pairCorr[i], rdf[i], idGasPairCorr[i]);
    }

    printf("\r\n");

    return lastError_;
}
