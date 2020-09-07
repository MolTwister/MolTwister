#include "CmdPotenEnergyMap.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../Utilities/DCDFile.h"
#include "../Tools/ProgressBar.h"

std::string CCmdPotenEnergyMap::getCmd()
{
    return "potenergymap";
}

std::vector<std::string> CCmdPotenEnergyMap::getCmdLineKeywords()
{
    return { "potenergymap", "coulomb" };
}

std::vector<std::string> CCmdPotenEnergyMap::getCmdHelpLines()
{
    return {
                "potenergymap <DCD filename> <frame from> <frame to> <atom IDs> <list of applied force fields> <Nx> <Ny> <cutting plane> <cutting plane pos>"
           };
}

std::string CCmdPotenEnergyMap::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates a 2D potential energy map by utilizing a test particle of charge +1, averaged over\r\n";
    text+= "\tframes from <frame from> to <frame to> in the DCD file specified by <DCD filename>. Only the\r\n";
    text+= "\tatoms listed in <atom IDs> (e.g., H, O, C7) yield a contribution to the map. The <atom IDs>\r\n";
    text+= "\tlist is to be comma separated without any spaces. The <list of applied forcefields> can be\r\n";
    text+= "\t* coulomb\r\n";
    text+= "\tThis list is also to be comma separated, without any spaces. <Nx> and <Ny> defines the number\r\n";
    text+= "\tof bins to include in the 2D map. The 2D map is calculated for a <cutting plane> which can be\r\n";
    text+= "\t* XY\r\n";
    text+= "\t* YZ\r\n";
    text+= "\t* ZX\r\n";
    text+= "\twhere the cutting plane position is <cutting plane pos>. If XY, then the cutting plane position\r\n";
    text+= "\tis along Z, if YZ, then it is along X and if ZX, then it is along Y.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. PlanePrinted : <cutting plane>(<0 if XY, 1 if YZ and 2 if ZX)\r\n";
    text+= "\t2. Selection : <atom ID 1> <atom ID 2> ... <atom ID n>\r\n";
    text+= "\t3. Size : (<Nx>, <Ny>)\r\n";
    text+= "\t4. DensityData={{<point x>, <point y>, <pot. E>}, {<point x>, <point y>, <pot. E>}, ..., {<point x>, <point y>, <pot. E>}};\r\n";
    text+= "\twhere <point x> and <point y> are points in the defined cutting plane.";

    return text;
}

std::string CCmdPotenEnergyMap::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<std::string> atomsToInclude;
    std::vector<std::string> ffToInclude;
    std::vector<int> atomIndicesToInclude;
    CDCDFile dcdFile;
    std::string text, cutString;
    bool includeCoulomb = false;
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


    // Find indices to loop over
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsToInclude = CASCIIUtility::getWords(text, ",");

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        text = state_->atoms_[i]->getID();
        for(int j=0; j<atomsToInclude.size(); j++)
        {
            if(text == atomsToInclude[j])
            {
                atomIndicesToInclude.emplace_back(i);
            }
        }
    }


    // Find which parts of force field to include
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    ffToInclude = CASCIIUtility::getWords(text, ",");

    for(int i=0; i<ffToInclude.size(); i++)
    {
        if(ffToInclude[i] == "coulomb") includeCoulomb = true;
        else
        {
            lastError_ = std::string("Error: do not recognize ") + ffToInclude[i] + std::string(" as valid part of force-field!");
            return lastError_;
        }
    }


    // Get specification of grid to sample
    C3DRect pbc = state_->view3D_->calcPBC();

    text = CASCIIUtility::getArg(arguments, arg++);
    int Na = atoi(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    int Nb = atoi(text.data());

    cutString = CASCIIUtility::getArg(arguments, arg++);
    int iCut; // 0=XY, 1=YZ, 2=ZX
    if(cutString == "XY") iCut = 0;
    else if(cutString == "YZ") iCut = 1;
    else if(cutString == "ZX") iCut = 2;
    else
    {
        lastError_ = std::string("Error: unkown cut ") + cutString + std::string("!");
        return lastError_;
    }


    if((Na <= 0) || (Nb <= 0))
    {
        lastError_ = std::string("Error cannot have zero or negative number of grid points (i.e. Na=") + std::to_string(Na) + std::string(", Nb=") + std::to_string(Nb) + std::string(" not possible)!");
        return lastError_;
    }

    double da=0.0, db=0.0, a_start=0.0, b_start=0.0;
    if(iCut == 0)
    {
        da = pbc.getWidthX() / double(Na);
        db = pbc.getWidthY() / double(Nb);
        a_start = pbc.rLow_.x_;
        b_start = pbc.rLow_.y_;
    }
    if(iCut == 1)
    {
        da = pbc.getWidthY() / double(Na);
        db = pbc.getWidthZ() / double(Nb);
        a_start = pbc.rLow_.y_;
        b_start = pbc.rLow_.z_;
    }
    if(iCut == 2)
    {
        da = pbc.getWidthZ() / double(Na);
        db = pbc.getWidthX() / double(Nb);
        a_start = pbc.rLow_.z_;
        b_start = pbc.rLow_.x_;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    double c_pos = atof(text.data());


    // Set up arrays to store results
    std::vector<double> points_a;
    std::vector<double> points_b;
    std::vector<double> values;
    points_a.resize(Na*Nb);
    points_b.resize(Na*Nb);
    values.resize(Na*Nb, 0.0);


    // Create map using minimum image convention
    const double K = 1389.35; // [(kJ/mol)AA]
    const double q_i = 1.0; // Partial charge of test particle
    int numFrames = dcdFile.getNumRecords();
    int count = 0;
    pb.beginProgress("Calculate potential energy map");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        C3DVector r_i;
        for(int ia=0; ia<Na; ia++)
        {
            double a_pos = a_start + double(ia)*da;
            if(iCut == 0) { r_i.x_ = a_pos; r_i.z_ = c_pos; }
            if(iCut == 1) { r_i.y_ = a_pos; r_i.x_ = c_pos; }
            if(iCut == 2) { r_i.z_ = a_pos; r_i.y_ = c_pos; }

            for(int ib=0; ib<Nb; ib++)
            {
                double b_pos = b_start + double(ib)*db;
                if(iCut == 0) { r_i.y_ = b_pos; }
                if(iCut == 1) { r_i.z_ = b_pos; }
                if(iCut == 2) { r_i.x_ = b_pos; }

                C3DVector displace;
                double Utot = 0.0;
                for(int imageX=-1; imageX<=1; imageX++)
                {
                    displace.x_ = double(imageX) * pbc.getWidthX();
                    for(int imageY=-1; imageY<=1; imageY++)
                    {
                        displace.y_ = double(imageY) * pbc.getWidthY();
                        for(int imageZ=-1; imageZ<=1; imageZ++)
                        {
                            displace.z_ = double(imageZ) * pbc.getWidthZ();
                            for(int j=0; j<atomIndicesToInclude.size(); j++)
                            {
                                C3DVector r_j = dcdFile.getCoordinate(atomIndicesToInclude[j]) + displace;

                                double r_ij = (r_i - r_j).norm();
                                if(r_ij <= 0.0) r_ij = 1E-9;

                                if(includeCoulomb)
                                {
                                    double q_j = state_->atoms_[atomIndicesToInclude[j]]->Q_;
                                    Utot+= K * (q_i * q_j) / r_ij;
                                }
                            }
                        }
                    }
                }

                int index = ia + ib*Na;
                if(t == frameFrom)
                {
                    points_a[index] = a_pos;
                    points_b[index] = b_pos;
                }
                values[index]+= Utot;
            }
        }

        count++;
        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();


    // Average over frames and print results
    if(count == 0)
    {
        lastError_ = std::string("Error: no frames found in range ") + std::to_string(frameFrom) + std::string(" to ") + std::to_string(frameTo) + std::string("!");
        return lastError_;
    }

    printf("\r\n");
    fprintf(stdOut_, "\tPlanePrinted : %s(%i)\r\n\tSelection :", cutString.data(), iCut);
    for(int i=0; i<atomsToInclude.size(); i++) fprintf(stdOut_, " %s", atomsToInclude[i].data());
    fprintf(stdOut_, "\r\n\tSize : (%i, %i)\r\n\tDensityData={", Na, Nb);
    for(int i=0; i<values.size(); i++)
    {
        values[i]/= double(count);

        if(i != 0) fprintf(stdOut_, ",");
        fprintf(stdOut_, "{%.6f, %.6f, %.10f}", points_a[i], points_b[i], values[i]);
    }
    fprintf(stdOut_, "};\r\n");

    return lastError_;
}
