#include "CmdFF.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdFF::getCmd()
{
    return "ff";
}

std::vector<std::string> CCmdFF::getCmdLineKeywords()
{
    return { "ff", "bondforceprofile", "angleforceprofile", "dihedralforceprofile", "nonbondforceprofile" };
}

std::vector<std::string> CCmdFF::getCmdHelpLines()
{
    return {
        "ff bondforceprofile <ff index> <profile start r> <profile end f> <num points in profile>",
        "ff angleforceprofile <ff index> <profile start theta> <profile end theta> <num points in profile>",
        "ff dihedralforceprofile <ff index> <profile start theta> <profile end theta> <num points in profile>",
        "ff nonbondforceprofile <ff index> <profile start r> <profile end f> <num points in profile>"
    };
}

std::string CCmdFF::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tThe molecular dynamics algorithm will first create a list of points based on the force-field\r\n";
    text+= "\tpotential expressions. This list of points will be interpolated during the simulations to\r\n";
    text+= "\tcalculate the forces that need to be appiled to each atom in the system.\r\n";
    text+= "\t\r\n";
    text+= "\tThis command will output the list of points that make up both the potential (in kJ/mol) and\r\n";
    text+= "\tthe force (in kJ/(mol*AA)). Depending on the type of potential we query (bond, angle, dihedral\r\n";
    text+= "\tor non-boned), the start and end points to output must be specified, together with the desired\r\n";
    text+= "\tnumber of points. The units for bond and non-bonded profiles are in Angstrom, while for angle\r\n";
    text+= "\tand dihedral profiles they are in degrees. The force-field index, <ff index>, can be found\r\n";
    text+= "\tthrough the 'list ff' command.";

    return text;
}

std::string CCmdFF::execute(std::vector<std::string> arguments)
{
    lastError_ = "";
    size_t arg = 0;

    std::string cmd =  CASCIIUtility::getArg(arguments, arg++).data();
    if(cmd == "bondforceprofile")
    {
        parseBondforceprofileCommand(arguments, arg);
    }
    else if(cmd == "angleforceprofile")
    {
        parseAngleforceprofileCommand(arguments, arg);
    }
    else if(cmd == "dihedralforceprofile")
    {
        parseDihedralforceprofileCommand(arguments, arg);
    }
    else if(cmd == "nonbondforceprofile")
    {
        parseNonbondforceprofileCommand(arguments, arg);
    }

    return lastError_;
}

void CCmdFF::parseBondforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg)
{
    int index = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    if((index < 0) || (index >= state_->mdFFBondList_.size()))
    {
        lastError_ = "Bond index " + std::to_string(index) + " does not exist!";
        return;
    }

    CMDFFBond* bond = state_->mdFFBondList_.get(index);
    if(!bond)
    {
        lastError_ = "Bond index " + std::to_string(index) + " was NULL!";
        return;
    }

    float rStart = (float)atof(CASCIIUtility::getArg(arguments, arg++).data());
    float rEnd = (float)atof(CASCIIUtility::getArg(arguments, arg++).data());
    int points = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    std::vector<std::pair<float, float>> forceProfile = bond->calc1DForceProfile(rStart, rEnd, points);
    std::vector<std::pair<float, float>> potentialProfile = bond->calc1DPotentialProfile(rStart, rEnd, points);

    fprintf(stdOut_, "\r\n\t%-15s%-15s%-15s\r\n", "r [Angstrom]", "Force", "Potential");
    for(size_t i=0; i<forceProfile.size(); i++)
    {
        if(i >= potentialProfile.size()) continue;
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", (double)forceProfile[i].first, (double)forceProfile[i].second, (double)potentialProfile[i].second);
    }
}

void CCmdFF::parseAngleforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg)
{
    int index = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    if((index < 0) || (index >= state_->mdFFAngleList_.size()))
    {
        lastError_ = "Angle index " + std::to_string(index) + " does not exist!";
        return;
    }

    CMDFFAngle* angle = state_->mdFFAngleList_.get(index);
    if(!angle)
    {
        lastError_ = "Angle index " + std::to_string(index) + " was NULL!";
        return;
    }

    float thetaStart = (float)atof(CASCIIUtility::getArg(arguments, arg++).data()) * float(M_PI / 180.0);
    float thetaEnd = (float)atof(CASCIIUtility::getArg(arguments, arg++).data()) * float(M_PI / 180.0);
    int points = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    std::vector<std::pair<float, float>> forceProfile = angle->calc1DForceProfile(thetaStart, thetaEnd, points);
    std::vector<std::pair<float, float>> potentialProfile = angle->calc1DPotentialProfile(thetaStart, thetaEnd, points);

    fprintf(stdOut_, "\r\n\t%-15s%-15s%-15s\r\n", "r [Angstrom]", "Force", "Potential");
    for(size_t i=0; i<forceProfile.size(); i++)
    {
        if(i >= potentialProfile.size()) continue;
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", (double)forceProfile[i].first * (180.0 / M_PI), (double)forceProfile[i].second, (double)potentialProfile[i].second);
    }
}

void CCmdFF::parseDihedralforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg)
{
    int index = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    if((index < 0) || (index >= state_->mdFFDihList_.size()))
    {
        lastError_ = "Dihedral index " + std::to_string(index) + " does not exist!";
        return;
    }

    CMDFFDih* dihedral = state_->mdFFDihList_.get(index);
    if(!dihedral)
    {
        lastError_ = "Dihedral index " + std::to_string(index) + " was NULL!";
        return;
    }

    float phiStart = (float)atof(CASCIIUtility::getArg(arguments, arg++).data()) * float(M_PI / 180.0);
    float phiEnd = (float)atof(CASCIIUtility::getArg(arguments, arg++).data()) * float(M_PI / 180.0);
    int points = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    std::vector<std::pair<float, float>> forceProfile = dihedral->calc1DForceProfile(phiStart, phiEnd, points);
    std::vector<std::pair<float, float>> potentialProfile = dihedral->calc1DPotentialProfile(phiStart, phiEnd, points);

    fprintf(stdOut_, "\r\n\t%-15s%-15s%-15s\r\n", "r [Angstrom]", "Force", "Potential");
    for(size_t i=0; i<forceProfile.size(); i++)
    {
        if(i >= potentialProfile.size()) continue;
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", (double)forceProfile[i].first * (180.0 / M_PI), (double)forceProfile[i].second, (double)potentialProfile[i].second);
    }
}

void CCmdFF::parseNonbondforceprofileCommand(const std::vector<std::string>& arguments, size_t& arg)
{
    int index = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    if((index < 0) || (index >= state_->mdFFNonBondedList_.size()))
    {
        lastError_ = "Non-bond index " + std::to_string(index) + " does not exist!";
        return;
    }

    CMDFFNonBonded* bond = state_->mdFFNonBondedList_.get(index);
    if(!bond)
    {
        lastError_ = "Non-bond index " + std::to_string(index) + " was NULL!";
        return;
    }

    float rStart = (float)atof(CASCIIUtility::getArg(arguments, arg++).data());
    float rEnd = (float)atof(CASCIIUtility::getArg(arguments, arg++).data());
    int points = atoi(CASCIIUtility::getArg(arguments, arg++).data());

    std::vector<std::pair<float, float>> forceProfile = bond->calc1DForceProfile(rStart, rEnd, points);
    std::vector<std::pair<float, float>> potentialProfile = bond->calc1DPotentialProfile(rStart, rEnd, points);

    fprintf(stdOut_, "\r\n\t%-15s%-15s%-15s\r\n", "r [Angstrom]", "Force", "Potential");
    for(size_t i=0; i<forceProfile.size(); i++)
    {
        if(i >= potentialProfile.size()) continue;
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", (double)forceProfile[i].first, (double)forceProfile[i].second, (double)potentialProfile[i].second);
    }
}
