#include "CmdFF.h"
#include "../../Utilities/ASCIIUtility.h"

std::string CCmdFF::getCmd()
{
    return "ff";
}

std::vector<std::string> CCmdFF::getCmdLineKeywords()
{
    return { "ff", "bondforceprofile", "angleforceprofile" };
}

std::vector<std::string> CCmdFF::getCmdHelpLines()
{
    return {
        "ff bondforceprofile <ff index> <profile start r> <profile end f> <num points in profile>",
        "ff angleforceprofile <ff index> <profile start theta> <profile end theta> <num points in profile>"
    };
}

std::string CCmdFF::getCmdFreetextHelp()
{
    std::string text;

    text+= "\t\r\n";

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
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", forceProfile[i].first, forceProfile[i].second, potentialProfile[i].second);
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
        fprintf(stdOut_, "\t% -15.4f% -15.4f% -15.4f\r\n", forceProfile[i].first * float(180.0 / M_PI), forceProfile[i].second, potentialProfile[i].second);
    }
}
