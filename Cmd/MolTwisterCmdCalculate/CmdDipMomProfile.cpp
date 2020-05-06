#include "CmdDipMomProfile.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/ProgressBar.h"

std::string CCmdDipMomProfile::getCmd()
{
    return "dipmomprofile";
}

std::vector<std::string> CCmdDipMomProfile::getCmdLineKeywords()
{
    return { "dipmomprofile", "chargedmol", "forcecomplete" };
}

std::vector<std::string> CCmdDipMomProfile::getCmdHelpLines()
{
    return {
                "dipmomprofile <DCD filename> <frame from> <frame to> <molecule resname> <x-dir.vec> <y-dir.vec> <z-dir.vec> <num. bins> <min. dist.> <max. dist.> [chargedmol] [forcecomplete <num. atoms in complete molecule>]"
           };
}

std::string CCmdDipMomProfile::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the dipole moment of the selected molecules from <molecule resname> (string\r\n";
    text+= "\twith resname). The dipolemoment is averaged based on all the defined molecules\r\n";
    text+= "\tin frames <frame from> to <frame to>. The DCD file, <DCD filename>, is used as input.\r\n";
    text+= "\t\r\n";
    text+= "\tThe dipole moment expression for neutral molecules are used as default. By sepcifying\r\n";
    text+= "\t'chargedmol', the dipole moment expression for charged molecules is employed. A profile\r\n";
    text+= "\tof the dipolemoment, for the defined molecule type, is created along the unit vector\r\n";
    text+= "\tcreated from the vector (<x-dir.vec>, <y-dir.vec>, <z-dir.vec>) from the distance <min.\r\n";
    text+= "\tdist.> to the distance <max. dist.>, where the distances are molecular center of mass\r\n";
    text+= "\t(COM) distances that are projected along the specified direction unit vector.The\r\n";
    text+= "\tprojected COM for each molecule where the dipole moment is calculated will be binned into\r\n";
    text+= "\t<num. bins> bins between the specified minimum and maximum distances. Each bin will contain\r\n";
    text+= "\tan average dipole moment from all molecules belonging to that bin. By default, the algorithm\r\n";
    text+= "\twill not check if a molecule is wrapped accross periodic boundaries, thus resulting in\r\n";
    text+= "\tincomplete molecules being processed. This can be fixed by applying the 'forcecomplete'\r\n";
    text+= "\tkeyword with the number of atoms to expect in a molecule, <num. atoms in complete molecule>.\r\n";
    text+= "\tMolecules that do not conform to the specified atom count will not be included in the\r\n";
    text+= "\tdipole moment profile.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\r\n";
    text+= "\t1. From frame: <first frame to include in average>\r\n";
    text+= "\t2. To frame: <last frame to include in average>\r\n";
    text+= "\t3. Resname: <resname for selected molecules>\r\n";
    text+= "\t4. Range: [<min. dist>, <max. dist>]\r\n";
    text+= "\t5. Num. bins: <num. bins>\r\n";
    text+= "\t6. Rc projection: <UnitVec(<x-dir.vec>, <y-dir.vec>, <z-dir.vec>)>\r\n";
    text+= "\t7. Distance Count P_x P_y P_z P_abs\r\n";
    text+= "\t8. <projected COM distance for bin (center of bin)> <num bin entries> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <absolute dipole moment>\r\n";
    text+= "\t9. <projected COM distance for bin (center of bin)> <num bin entries> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <absolute dipole moment>\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\tN+7. <projected COM distance for bin (center of bin)> <num bin entries> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <absolute dipole moment>\r\n";
    text+= "\twhere N is the number of bins.";

    return text;
}

std::string CCmdDipMomProfile::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<C3DVector> profile;
    std::vector<double> absDipMom;
    std::vector<int> countArray;
    std::vector<std::vector<int>> molecules;
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


    // Obtain list of all molecules to include in dipole moment calculations
    text = CASCIIUtility::getArg(arguments, arg++);
    std::string resnameString = text;

    CMolTwisterStateTools stateTools(state_, stdOut_);
    int numMolecules = stateTools.getNumMolecules();
    for(int i=0; i<numMolecules; i++)
    {
        std::vector<int> aAtomIndices;
        stateTools.getAtomsOfMolecule(i, aAtomIndices);
        if(aAtomIndices.size() <= 0) continue;

        text = state_->atoms_[aAtomIndices[0]]->resname_;
        if(text == resnameString)
        {
            molecules.emplace_back(aAtomIndices);
        }
    }


    // Obtain vector to project molecular center of mass (Rc) along
    C3DVector vecProjection;
    text = CASCIIUtility::getArg(arguments, arg++);
    vecProjection.x_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecProjection.y_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vecProjection.z_ = atof(text.data());
    C3DVector vecProjectionUnit = vecProjection.unit();


    // Obtain number of profile bins and profile boundaries
    text = CASCIIUtility::getArg(arguments, arg++);

    int size = atoi(text.data());
    if(size <= 0)
    {
        lastError_ = "Error: cannot have zero or fewer bins in profile!";
        return lastError_;
    }

    profile.resize(size, C3DVector(0.0, 0.0, 0.0));
    countArray.resize(size, 0);
    absDipMom.resize(size, 0.0);

    text = CASCIIUtility::getArg(arguments, arg++);
    double minDist = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double maxDist = atof(text.data());

    if(fabs(minDist - maxDist) < double(FLT_MIN))
    {
        lastError_ = "Error: cannot have equal minimum and maximum distances in profile!";
        return lastError_;
    }


    // Are we going to use charge neutral formulation: sum_i q_i r_i or
    // are we using charged molecules: sum_i q_i (r_i - R_c)?
    bool chargeNeutralFormulation = true;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "chargedmol")
        chargeNeutralFormulation = false;
    else arg--;


    // Force to only consider complete molecules, specified
    // by the number of atoms within the molecule (can be
    // used instead of performing atomic unwrap of a system)
    int numAtomsInMol = -1;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "forcecomplete")
    {
        text = CASCIIUtility::getArg(arguments, arg++);

        numAtomsInMol = atoi(text.data());
    }
    else arg--;


    // Go through all frames and build profile
    int numFrames = dcdFile.getNumRecords();
    pb.beginProgress("Calculate dipole moment profile");
    for(int t=frameFrom; t<frameTo; t++)
    {
        if((t < 0) || (t >= numFrames))
        {
            lastError_ = std::string("Error: frame ") + std::to_string(t) + std::string(" does not exist (num. frames = ") + std::to_string(numFrames) + std::string(")");
            return lastError_;
        }

        dcdFile.gotoRecord(t);

        C3DVector P, Rc;
        CDCDTools dcdTools(state_, stdOut_);
        for(int i=0; i<molecules.size(); i++)
        {
            if((numAtomsInMol != -1) && (molecules[i].size() != numAtomsInMol))
                continue;

            P = dcdTools.getMoleculeDipoleMoment(molecules[i], &dcdFile, Rc, chargeNeutralFormulation);
            double projDist = Rc * vecProjectionUnit;

            int n = int((projDist - minDist) / (maxDist - minDist) * double(size));
            if((n >= 0) && (n < size))
            {
                profile[n]+= P;
                absDipMom[n]+= P.norm();
                countArray[n]++;
            }
        }

        pb.updateProgress(t-frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();

    for(int i=0; i<size; i++)
    {
        double norm = (countArray[i] == 0) ? 1.0 : (1.0 / double(countArray[i]));
        profile[i]*= norm;
        absDipMom[i]*= norm;
    }


    // Print results
    printf("\r\n");
    fprintf(stdOut_, "From frame: %i\r\n", frameFrom);
    fprintf(stdOut_, "To frame: %i\r\n", frameTo);
    fprintf(stdOut_, "Resname: %s\r\n", resnameString.data());
    fprintf(stdOut_, "Range: [%.4f, %.4f]\r\n", minDist, maxDist);
    fprintf(stdOut_, "Num. bins: %i\r\n\r\n", size);
    fprintf(stdOut_, "Rc projection: (%g, %g, %g)\r\n\r\n", vecProjectionUnit.x_, vecProjectionUnit.y_, vecProjectionUnit.z_);
    fprintf(stdOut_, "%-15s%-15s%-15s%-15s%-15s%-15s\r\n", "Distance", "Count", "P_x", "P_y", "P_z", "P_abs");
    fprintf(stdOut_, "---------------------------------------------\r\n");

    for(int i=0; i<profile.size(); i++)
    {
        double dDist = 0.5 * double(2*i+1) * (maxDist - minDist) / double(size) + minDist;
        fprintf(stdOut_, "%-15.4f%-15i%-15.6f%-15.6f%-15.6f%-15.6f\r\n", dDist, countArray[i], profile[i].x_, profile[i].y_, profile[i].z_, absDipMom[i]);
    }

    printf("\r\n");

    return lastError_;
}
