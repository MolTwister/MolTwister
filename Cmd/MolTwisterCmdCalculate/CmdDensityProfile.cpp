#include "CmdDensityProfile.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/ProgressBar.h"
#include "../Tools/DCDTools.h"

std::string CCmdDensityProfile::getCmd()
{
    return "densityprofile";
}

std::vector<std::string> CCmdDensityProfile::getCmdLineKeywords()
{
    return { "densityprofile", "within", "simbox", "cylinder" };
}

std::vector<std::string> CCmdDensityProfile::getCmdHelpLines()
{
    return {
                "densityprofile <DCD filename> <frame from> <frame to> <atoms (comma sep, no space)> <vec. from> <vec. to> <num bins> [within simbox <vec. low> <vec. high>, within cylinder <radius>]"
           };
}

std::string CCmdDensityProfile::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the density profile from an atomic selection, sepcified by the list of atom IDs,\r\n";
    text+= "\t<atoms>, s.a., O, H, C5, etc (no cammas). The atomic positions are taken from the specified\r\n";
    text+= "\tDCD file, from <frame from> to <frame to>. The density profile is calculated along the vector\r\n";
    text+= "\tV = <vec. to> - <vec. from>, where <vec. from> and <vec. to>, both are specified on the form\r\n";
    text+= "\t<x> <y> <z>. The density profile is divided into <num bins> between these two points. Optionally,\r\n";
    text+= "\tthe 'within' keyword can be used to specify a cutoff for which atoms are not included into the\r\n";
    text+= "\tdensity profile. If simbox is used, then atoms outside the lower corner <vec. low> and the upper\r\n";
    text+= "\tcorner <vec. high> (both of the form <x> <y> <z>), are not counted. If cylinder is used, then all\r\n";
    text+= "\tatoms outside the cylinder of radius <radius> with its core along V are ignored. If within is not\r\n";
    text+= "\tspecified, then all atoms are counted. The output is specified below.\r\n";
    text+= "\r\n";
    text+= "\tNote that there are three types of calculated densities from this command, these are\r\n";
    text+= "\t* Particle density\r\n";
    text+= "\t* Mass density\r\n";
    text+= "\t* Charge density\r\n";
    text+= "\r\n";
    text+= "\tThe output is as follows, starting with a single header line,\r\n";
    text+= "\t1. abs_x abs_y abs_z rel_x rel_y rel_z #/frm m/frm Q/frm\r\n";
    text+= "\t2. <x abs> <y <bs> <z abs> <x rel> <y rel> <z rel> <particle density> <mass density> <charge density>\r\n";
    text+= "\t3. <x abs> <y <bs> <z abs> <x rel> <y rel> <z rel> <particle density> <mass density> <charge density>\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\tN+1. <x abs> <y <bs> <z abs> <x rel> <y rel> <z rel> <particle density> <mass density> <charge density>\r\n";
    text+= "\twhere N is the number of bins. Relative positions (rel) are relative to <vec. from>, while absolute\r\n";
    text+= "\tpositions (abs) are relative to the simulation box. Mass density is in [g/mol], while charge density\r\n";
    text+= "\tis in units of partial charges.";

    return text;
}

std::string CCmdDensityProfile::execute(std::vector<std::string> arguments)
{
    /////////////////////////////////////////////////////////////
    // Let $\mathbf{v}_f$ and $\mathbf{v}_t$ be the from and to vectors, respectively, that make up the vector
    // $$\mathbf{k} = \mathbf{v}_t - \mathbf{v}_f,$$
    // which represent the line through 3D-space in which we want to estimate $N$ evenly spaced densities.
    //
    // Let $\mathbf{r}$ be the position of some particle that we want to project to $\mathbf{k}$ and let $\mathbf{R}$ be the particle position relative to $\mathbf{v}_f$. Then,
    // $$\mathbf{R} = \mathbf{r} - \mathbf{v}_f.$$
    // The projection, $p_k$, of the particle at $\mathbf{r}$ (i.e. the length from $\mathbf{v}_f$ to the point at $\mathbf{k}$ which makes a perpendicular line to $\mathbf{r}$) can be found be
    // $$\cos \theta = \frac{\mathbf{k}\cdot \mathbf{R}}{\vert \mathbf{k} \vert \vert \mathbf{R} \vert} = \frac{p_k}{\vert \mathbf{R} \vert}.$$
    // Thus, $\mathbf{p}_k = \mathbf{R}\cdot \hat{mathbf{k}}$, where $\hat{mathbf{k}}$ is a unit vector.
    /////////////////////////////////////////////////////////////

    lastError_ = "";

    size_t arg = 0;
    CProgressBar pb;
    std::vector<C3DVector> absPosAlong_k;
    std::vector<C3DVector> relPosAlong_k;
    std::vector<std::string> atomsToInclude;
    std::vector<double> atomicMasses;
    std::vector<double> atomicCharges;
    std::vector<double> massDensitiesAlong_k;
    std::vector<double> chargeDensitiesAlong_k;
    std::vector<int> atomIndicesToInclude;
    std::vector<int> numDensitiesAlong_k;
    C3DVector vf;
    C3DVector vt;
    CDCDFile dcdFile;
    std::string text;
    int numFrames, frameFrom, frameTo, count;


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


    // Find indices to loop over and their masses
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsToInclude = CASCIIUtility::getWords(text.data(), ",");

    for(int i=0; i<state_->atoms_.size(); i++)
    {
        text = state_->atoms_[i]->getID();
        for(int j=0; j<atomsToInclude.size(); j++)
        {
            if(text == atomsToInclude[j])
            {
                atomIndicesToInclude.emplace_back(i);
                atomicMasses.emplace_back(state_->atoms_[i]->m_);
                atomicCharges.emplace_back(state_->atoms_[i]->Q_);
            }
        }
    }


    // Calculate vector k to bin particles along
    text = CASCIIUtility::getArg(arguments, arg++);
    vf.x_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vf.y_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vf.z_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    vt.x_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vt.y_ = atof(text.data());
    text = CASCIIUtility::getArg(arguments, arg++);
    vt.z_ = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    numDensitiesAlong_k.resize(atoi(text.data()), 0);
    massDensitiesAlong_k.resize(atoi(text.data()), 0.0);
    chargeDensitiesAlong_k.resize(atoi(text.data()), 0.0);

    int withinType = 0;
    double withinCylRadius = 0.0;
    C3DRect withinSimBox;
    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "within")
    {
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text == "simbox")
        {
            withinType = 1;

            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rLow_.x_ = atof(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rLow_.y_ = atof(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rLow_.z_ = atof(text.data());

            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rHigh_.x_ = atof(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rHigh_.y_ = atof(text.data());
            text = CASCIIUtility::getArg(arguments, arg++);
            withinSimBox.rHigh_.z_ = atof(text.data());
        }
        else if(text == "cylinder")
        {
            withinType = 2;

            text = CASCIIUtility::getArg(arguments, arg++);
            withinCylRadius = atof(text.data());
        }
        else
        {
            lastError_ = "Error: expected argument to be 'simbox' or 'cylinder' after 'within'!";
            return lastError_;
        }
    }
    else arg--;

    int n;
    double pk, len_k, alpha;
    C3DVector r, R, pos;
    C3DVector vk = vt - vf;
    C3DVector vku = vk.unit();


    // Initialize vectors to contain the profile
    for(int i=0; i<numDensitiesAlong_k.size(); i++)
    {
        alpha = 0.5*( (2.0*double(i) + 1.0) / double(numDensitiesAlong_k.size()) );
        pos = vk*alpha;

        numDensitiesAlong_k[i] = 0;
        relPosAlong_k.emplace_back(pos);
        absPosAlong_k.emplace_back(vf + pos);
    }


    // Go through each frame and bin particles along vector k
    len_k = vk.norm();
    numFrames = dcdFile.getNumRecords();
    count = 0;
    pb.beginProgress("Calculating densities");
    for(int i=frameFrom; i<frameTo; i++)
    {
        if((i < 0) || (i >= numFrames)) continue;

        dcdFile.gotoRecord(i);

        for(int j=0; j<atomIndicesToInclude.size(); j++)
        {
            r = dcdFile.getCoordinate(atomIndicesToInclude[j]);
            R = r - vf;
            pk = R*vku;

            bool valid = true;
            if(withinType == 1)
            {
                if(r.x_ < withinSimBox.rLow_.x_) valid = false;
                if(r.y_ < withinSimBox.rLow_.y_) valid = false;
                if(r.z_ < withinSimBox.rLow_.z_) valid = false;

                if(r.x_ > withinSimBox.rHigh_.x_) valid = false;
                if(r.y_ > withinSimBox.rHigh_.y_) valid = false;
                if(r.z_ > withinSimBox.rHigh_.z_) valid = false;
            }
            if(withinType == 2)
            {
                C3DVector p_perp = R - vku*pk;
                double dp_perp = p_perp.norm();
                if(dp_perp > withinCylRadius) valid = false;
            }

            n = int((pk / len_k)*double(numDensitiesAlong_k.size()));
            if(valid && (n >= 0) && (n < numDensitiesAlong_k.size()))
            {
                numDensitiesAlong_k[n]++;
                massDensitiesAlong_k[n]+= atomicMasses[j]; // g/mol
                chargeDensitiesAlong_k[n]+= atomicCharges[j]; // Partial charges
            }
        }

        count++;

        pb.updateProgress(i - frameFrom, frameTo - frameFrom);
    }
    pb.endProgress();

    if(count == 0)
    {
        lastError_ = std::string("Error: no frames in range [") + std::to_string(frameFrom) + std::string(",") + std::to_string(frameTo) + std::string(")!");
        return lastError_;
    }


    // Print the results
    fprintf(stdOut_, "\t%-10s%-10s%-10s%-10s%-10s%-10s%-15s%-15s%-15s\r\n", "abs_x", "abs_y", "abs_z", "rel_x", "rel_y", "rel_z", "#/frm", "m/frm", "Q/frm");
    for(int i=0; i<numDensitiesAlong_k.size(); i++)
    {
        double avg_cnt = double(numDensitiesAlong_k[i]) / double(count);
        double avg_m = massDensitiesAlong_k[i] / double(count);
        double avg_Q = chargeDensitiesAlong_k[i] / double(count);
        fprintf(stdOut_, "\t% -10.4f% -10.4f% -10.4f% -10.4f% -10.4f% -10.4f%-15.6f%-15.6f%-15.6f\r\n",
                absPosAlong_k[i].x_, absPosAlong_k[i].y_, absPosAlong_k[i].z_,
                relPosAlong_k[i].x_, relPosAlong_k[i].y_, relPosAlong_k[i].z_, avg_cnt, avg_m, avg_Q);
    }

    return lastError_;
}
