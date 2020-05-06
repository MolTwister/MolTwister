#include "CmdCoulombEnergy.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/MolTwisterStateTools.h"
#include "../Tools/ParsingTools.h"
#include "../Tools/MolecularTools.h"

std::string CCmdCoulombEnergy::getCmd()
{
    return "coulombenergy";
}

std::vector<std::string> CCmdCoulombEnergy::getCmdLineKeywords()
{
    return { "coulombenergy", "single", "dihedralrot", "anglerot", "bondstretch", "allframes", "rot", "stretch", "id", "var", "1to4coeffs" };
}

std::vector<std::string> CCmdCoulombEnergy::getCmdHelpLines()
{
    return {
                "coulombenergy single [1to4coeffs <coeff 1> <coeff 2> <coeff 3> <coeff 4>]",
                "coulombenergy dihedralrot <dihedral> rot <start angle> <end angle> <angular step size> [1to4coeffs <coeff 1> <coeff 2> <coeff 3> <coeff 4>]",
                "coulombenergy anglerot <angle> rot <start angle> <end angle> <angular step size> [1to4coeffs <coeff 1> <coeff 2> <coeff 3> <coeff 4>]",
                "coulombenergy bondstretch <bond> stretch <start dist> <end dist> <dist step size> [1to4coeffs <coeff 1> <coeff 2> <coeff 3> <coeff 4>]",
                "coulombenergy allframes [1to4coeffs <coeff 1> <coeff 2> <coeff 3> <coeff 4>]"
           };
}

std::string CCmdCoulombEnergy::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tMeasures the Coulomb energy using one of the following choices\r\n";
    text+= "\t* 'single':      sum Coulomb energies between all atoms, N, by (i=0..N, j=i+1..N)\r\n";
    text+= "\t                 with output 'Etot = <energy> kJ/mol'.\r\n";
    text+= "\t* 'dihedralrot': performs a 'single' for dihedral angles, specified by the\r\n";
    text+= "\t                 <dihedral>, <start angle>, <end angle>, <angular step size>\r\n";
    text+= "\t                 parameters. The output is a two column space sepatated list\r\n";
    text+= "\t                 with header 'Angle [deg] Etot [kJ/mol]'.\r\n";
    text+= "\t* 'anglerot':    performs a 'single' for molecular angles, specified by the\r\n";
    text+= "\t                 <angle>, <start angle>, <end angle>, <angular step size>\r\n";
    text+= "\t                 parameters. The output is a two column space sepatated list\r\n";
    text+= "\t                 with header 'Angle [deg] Etot [kJ/mol]'.\r\n";
    text+= "\t* 'bondstretch': performs a 'single' for bond stretching, specified by the\r\n";
    text+= "\t                 <bond>, <start dist>, <end dist>, <dist step size>\r\n";
    text+= "\t                 parameters. The output is a two column space sepatated list\r\n";
    text+= "\t                 with header 'Dist [Å] Etot [kJ/mol]'.\r\n";
    text+= "\t* 'allframes':   performs a 'single' for all the loaded frames and produces\r\n";
    text+= "\t                 a two column output of the form 'Frame Etot [kJ/mol].'\r\n";
    text+= "\r\n";
    text+= "\tThe <dihedral> parameter specifies the dihedral to rotate and can be formatted\r\n";
    text+= "\tas shown below.\r\n";
    text+= "\t* id <atom index 1> <atom index 2> <atom index 3> <atom index 4>\r\n";
    text+= "\t* var <variable name>\r\n";
    text+= "\r\n";
    text+= "\tThe <angle> parameter specifies the angle to rotate and can be formatted as\r\n";
    text+= "\tshown below.\r\n";
    text+= "\t* id <atom index 1> <atom index 2> <atom index 3>\r\n";
    text+= "\t* var <variable name>\r\n";
    text+= "\r\n";
    text+= "\tThe <bond> parameter specifies the bond to stretch and can be formatted as\r\n";
    text+= "\tshown below.\r\n";
    text+= "\t* id <atom index 1> <atom index 2>\r\n";
    text+= "\t* var <variable name>\r\n";
    text+= "\t\r\n";
    text+= "\tThe '1to4coeffs' keyword can be applied to assign different weights to the\r\n";
    text+= "\tCoulomb energy calculations that are separated with 1, 2, 3 and 4 bonds,\r\n";
    text+= "\trespectively, by assigning values to <coeff 1> through <coeff 4>.";

    return text;
}

std::string CCmdCoulombEnergy::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    CMolTwisterStateTools stateTools(state_, stdOut_);
    double a1to4BondSepCoeffs[4];
    std::string text;
    int pathID = -1;

    text = CASCIIUtility::getArg(arguments, arg++);
    if(text == "single")
    {
        pathID = 1;
    }
    else if(text == "dihedralrot")
    {
        pathID = 2;
    }
    else if(text == "anglerot")
    {
        pathID = 3;
    }
    else if(text == "bondstretch")
    {
        pathID = 4;
    }
    else if(text == "allframes")
    {
        pathID = 5;
    }
    else
    {
        lastError_ = "Syntax Error: Third argument should indicate for example single point, dihedral rotation, etc.!";
        return lastError_;
    }

    if(pathID == 1)
    {
        double E;

        if(!CParsingTools(state_, stdOut_).retrieve1to4BondSepCoeffs(arguments, arg, a1to4BondSepCoeffs))
            E = stateTools.measureTotCoulombEnergy();
        else
            E = stateTools.measureTotCoulombEnergy(a1to4BondSepCoeffs);

        fprintf(stdOut_, "\r\nEtot = %g kJ/mol\r\n", E);
    }
    else if(pathID == 2)
    {
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        CAtom* atom3Ptr = nullptr;
        CAtom* atom4Ptr = nullptr;
        double E;
        bool special1to4;
        int atomIndices[4];

        std::pair<bool, std::string> retVal = CParsingTools(state_, stdOut_).retrieveDihedralAtoms(arguments, arg, &atom1Ptr, &atom2Ptr, &atom3Ptr, &atom4Ptr, atomIndices);
        if(retVal.first)
        {
            double startAngle, endAngle, step;

            text = CASCIIUtility::getArg(arguments, arg++);
            if(text == "rot")
            {
                text = CASCIIUtility::getArg(arguments, arg++);
                startAngle = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                endAngle = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                step = atof(text.data());

                special1to4 = CParsingTools(state_, stdOut_).retrieve1to4BondSepCoeffs(arguments, arg, a1to4BondSepCoeffs);

                fprintf(stdOut_, "\r\n\t%-15s%-15s\r\n\t------------------------------------\r\n", "Angle [Deg]", "Etot [kJ/mol]");
                if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                {
                    for(double angle=startAngle; angle<=endAngle; angle+=step)
                    {
                        CMolecularTools::modDihedralTo(atom1Ptr, atom2Ptr, atom3Ptr, atom4Ptr, angle * M_PI/180.0, state_->getCurrFrameIndex());
                        if(!special1to4)  E = stateTools.measureTotCoulombEnergy();
                        else              E = stateTools.measureTotCoulombEnergy(a1to4BondSepCoeffs);

                        fprintf(stdOut_, "\t% -15g% -15g\r\n", angle, E);
                    }
                    state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                }
            }
            else
            {
                lastError_ = "Syntax Error: Fifth argument should be 'rot'!";
                return lastError_;
            }
        }
        else
        {
            lastError_ = retVal.second;
            return lastError_;
        }
    }
    else if(pathID == 3)
    {
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        CAtom* atom3Ptr = nullptr;
        double E;
        bool special1to4;
        int atomIndices[3];

        std::pair<bool, std::string> retVal = CParsingTools(state_, stdOut_).retrieveAngleAtoms(arguments, arg, &atom1Ptr, &atom2Ptr, &atom3Ptr, atomIndices);
        if(retVal.first)
        {
            double startAngle, endAngle, step;

            text = CASCIIUtility::getArg(arguments, arg++);
            if(text == "rot")
            {
                text = CASCIIUtility::getArg(arguments, arg++);
                startAngle = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                endAngle = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                step = atof(text.data());

                special1to4 = CParsingTools(state_, stdOut_).retrieve1to4BondSepCoeffs(arguments, arg, a1to4BondSepCoeffs);

                fprintf(stdOut_, "\r\n\t%-15s%-15s\r\n\t------------------------------------\r\n", "Angle [Deg]", "Etot [kJ/mol]");
                if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                {
                    for(double angle=startAngle; angle<=endAngle; angle+=step)
                    {
                        CMolecularTools::modAngleTo(atom1Ptr, atom2Ptr, atom3Ptr, angle * M_PI/180.0, state_->getCurrFrameIndex());
                        if(!special1to4)   E = stateTools.measureTotCoulombEnergy();
                        else               E = stateTools.measureTotCoulombEnergy(a1to4BondSepCoeffs);

                        fprintf(stdOut_, "\t% -15g% -15g\r\n", angle, E);
                    }
                    state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                }
            }
            else
            {
                lastError_ = "Syntax Error: Fifth argument should be 'rot'!";
                return lastError_;
            }
        }
        else
        {
            lastError_ = retVal.second;
            return lastError_;
        }
    }
    else if(pathID == 4)
    {
        CAtom* atom1Ptr = nullptr;
        CAtom* atom2Ptr = nullptr;
        double E;
        bool special1to4;
        int atomIndices[2];

        std::pair<bool, std::string> retVal = CParsingTools(state_, stdOut_).retrieveBondAtoms(arguments, arg, &atom1Ptr, &atom2Ptr, atomIndices);
        if(retVal.first)
        {
            double startDist, endDist, step;

            text = CASCIIUtility::getArg(arguments, arg++);
            if(text == "stretch")
            {
                text = CASCIIUtility::getArg(arguments, arg++);
                startDist = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                endDist = atof(text.data());
                text = CASCIIUtility::getArg(arguments, arg++);
                step = atof(text.data());

                special1to4 = CParsingTools(state_, stdOut_).retrieve1to4BondSepCoeffs(arguments, arg, a1to4BondSepCoeffs);

                fprintf(stdOut_, "\r\n\t%-15s%-15s\r\n\t------------------------------------\r\n", "Dist [Å]", "Etot [kJ/mol]");
                if(state_->saveCoordinates(state_->getCurrFrameIndex()))
                {
                    for(double dist=startDist; dist<=endDist; dist+=step)
                    {
                        CMolecularTools::modBondLengthTo(atom1Ptr, atom2Ptr, dist, state_->getCurrFrameIndex());
                        if(!special1to4)   E = stateTools.measureTotCoulombEnergy();
                        else               E = stateTools.measureTotCoulombEnergy(a1to4BondSepCoeffs);

                        fprintf(stdOut_, "\t% -15g% -15g\r\n", dist, E);
                    }
                    state_->retrieveSavedCoordinates(state_->getCurrFrameIndex());
                }
            }
            else
            {
                lastError_ = "Syntax Error: Fifth argument should be 'rot'!";
                return lastError_;
            }
        }
        else
        {
            lastError_ = retVal.second;
            return lastError_;
        }
    }
    else if(pathID == 5)
    {
        double E;
        bool special1to4;

        if(state_->atoms_.size() <= 0)
        {
            lastError_ = "No atoms found!";
            return lastError_;
        }

        special1to4 = CParsingTools(state_, stdOut_).retrieve1to4BondSepCoeffs(arguments, arg, a1to4BondSepCoeffs);

        fprintf(stdOut_, "\r\n\t%-15s%-15s\r\n\t------------------------------------\r\n", "Frame", "Etot [kJ/mol]");
        for(int i=0; i<state_->atoms_[0]->r_.size(); i++)
        {
            if(!special1to4)   E = stateTools.measureTotCoulombEnergy(nullptr, i);
            else               E = stateTools.measureTotCoulombEnergy(a1to4BondSepCoeffs, i);

            fprintf(stdOut_, "\t% -15i% -15g\r\n", i, E);
        }
    }
    else
    {
        lastError_ = "Error: Unrecognized parameters!";
    }

    return lastError_;
}
