#include "CmdForceBetween.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../../MDFF/MolTwisterMDFFCoulomb.h"

std::string CCmdForceBetween::getCmd()
{
    return "forcebetween";
}

std::vector<std::string> CCmdForceBetween::getCmdLineKeywords()
{
    return { "forcebetween", "nonbonded", "coulomb", "bond", "angle", "dihedral", "numerical" };
}

std::vector<std::string> CCmdForceBetween::getCmdHelpLines()
{
    return {
                "forcebetween nonbonded <comma sep. list of atom IDs> <FF index (zero based)> [numerical]",
                "forcebetween coulomb <comma sep. list of atom IDs>",
                "forcebetween bond <comma sep. list of atom IDs> <FF index (zero based)> [numerical]",
                "forcebetween angle <comma sep. list of atom IDs> <FF index (zero based)> [numerical]",
                "forcebetween dihedral <comma sep. list of atom IDs> <FF index (zero based)> [numerical]"
           };
}

std::string CCmdForceBetween::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculates the nonbonded, Coulomb or bond force (in kJ/(mol Angstrom)) betwee two atoms, based on\r\n";
    text+= "\tthe configured force fields, or the angular force between three atoms or the dihedral force\r\n";
    text+= "\tbetween four atoms. The forces are calculated based on the atoms shown in the current frame.\r\n";
    text+= "\tThe indices of the various force fields that have been created or loaded is available through\r\n";
    text+= "\tthe 'list' command. An index into these lists, for the appropriate force field type, is specified\r\n";
    text+= "\tthrough the <FF index (zero based)> parameter. The force is calculated between the atoms given\r\n";
    text+= "\tin <comma sep. list of atom IDs> (e.g., H, O, C7, where the list must be enetered withour spece).\r\n";
    text+= "\tBy sepcifying 'numerical', the forces are calculated using numerical differentiation.";

    return text;
}

std::string CCmdForceBetween::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text, forceTypeString;
    std::vector<std::string> atomsToInclude;


    // Get force type
    forceTypeString = CASCIIUtility::getArg(arguments, arg++);

    // Find atom indices to calculate forces between
    text = CASCIIUtility::getArg(arguments, arg++);

    CASCIIUtility::removeWhiteSpace(text);
    atomsToInclude = CASCIIUtility::getWords(text, ",");

    // Calculate forces
    if(forceTypeString == "nonbonded")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for non-bonded force calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFNonBondedList_.size()))
        {
            lastError_ = std::string("Error: non-bonded force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        bool numerical = false;
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "numerical") numerical = true;

        C3DVector f1, f2;
        if(!state_->mdFFNonBondedList_.get(ffIndex)->isPotentialCalcAvailable() ||
           !state_->mdFFNonBondedList_.get(ffIndex)->isForceCalcAvailable())
        {
            lastError_ = std::string("Error: non-bonded force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of force calculations!");
            return lastError_;
        }
        if(numerical)
        {
            if(!state_->mdFFNonBondedList_.get(ffIndex)->calcForcesNumerically(r1, r2, f1, f2))
            {
                printf("Warning: numerical derivatives failed to converge!\r\n");
            }
        }
        else
        {
            state_->mdFFNonBondedList_.get(ffIndex)->calcForces(r1, r2, f1, f2);
        }

        printf("\r\n");
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f1 = (%g,%g,%g), |f1| = %g\r\n", state_->mdFFNonBondedList_.get(ffIndex)->getFFType().data(), atomName1.data(), index1, r1.x_, r1.y_, r1.z_, f1.x_, f1.y_, f1.z_, f1.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f2 = (%g,%g,%g), |f2| = %g\r\n", state_->mdFFNonBondedList_.get(ffIndex)->getFFType().data(), atomName2.data(), index2, r2.x_, r2.y_, r2.z_, f2.x_, f2.y_, f2.z_, f2.norm());
        if(numerical) fprintf(stdOut_, "\r\n\tNote: calculations were done numerically!\r\n");
    }
    else if(forceTypeString == "coulomb")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for Coulomb force calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        double q1 = state_->atoms_[index1]->Q_;
        double q2 = state_->atoms_[index2]->Q_;
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        bool numerical = false;
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "numerical") numerical = true;

        C3DVector f1, f2;
        CMDFFCoulomb Coulomb(state_);
        if(numerical)
        {
            double distR12;
            if(!Coulomb.calcForceBetweenNumerically(r1, r2, q1, q2, f1, f2, distR12))
            {
                printf("Warning: numerical derivatives failed to converge!\r\n");
            }
        }
        else
        {
            Coulomb.calcForceBetween(r1, r2, q1, q2, f1, f2);
        }

        printf("\r\n");
        fprintf(stdOut_, "\tForce (Coulomb) on %s (index %i) at pos (%g,%g,%g): f1 = (%g,%g,%g), |f1| = %g\r\n", atomName1.data(), index1, r1.x_, r1.y_, r1.z_, f1.x_, f1.y_, f1.z_, f1.norm());
        fprintf(stdOut_, "\tForce (Coulomb) on %s (index %i) at pos (%g,%g,%g): f2 = (%g,%g,%g), |f2| = %g\r\n", atomName2.data(), index2, r2.x_, r2.y_, r2.z_, f2.x_, f2.y_, f2.z_, f2.norm());
        if(numerical) fprintf(stdOut_, "\r\n\tNote: calculations were done numerically!\r\n");
    }
    else if(forceTypeString == "bond")
    {
        if(atomsToInclude.size() != 2)
        {
            lastError_ = "Error: only two atoms can be specified for bond force calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFBondList_.size()))
        {
            lastError_ = std::string("Error: bond force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();

        bool numerical = false;
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "numerical") numerical = true;

        C3DVector f1, f2;
        if(!state_->mdFFBondList_.get(ffIndex)->isPotentialCalcAvailable() ||
           !state_->mdFFBondList_.get(ffIndex)->isForceCalcAvailable())
        {
            lastError_ = std::string("Error: bond force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of force calculations!");
            return lastError_;
        }
        if(numerical)
        {
            if(!state_->mdFFBondList_.get(ffIndex)->calcForcesNumerically(r1, r2, f1, f2))
            {
                printf("Warning: numerical derivatives failed to converge!\r\n");
            }
        }
        else
        {
            state_->mdFFBondList_.get(ffIndex)->calcForces(r1, r2, f1, f2);
        }

        printf("\r\n");
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f1 = (%g,%g,%g), |f1| = %g\r\n", state_->mdFFBondList_.get(ffIndex)->getFFType().data(), atomName1.data(), index1, r1.x_, r1.y_, r1.z_, f1.x_, f1.y_, f1.z_, f1.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f2 = (%g,%g,%g), |f2| = %g\r\n", state_->mdFFBondList_.get(ffIndex)->getFFType().data(), atomName2.data(), index2, r2.x_, r2.y_, r2.z_, f2.x_, f2.y_, f2.z_, f2.norm());
        if(numerical) fprintf(stdOut_, "\r\n\tNote: calculations were done numerically!\r\n");
    }
    else if(forceTypeString == "angle")
    {
        if(atomsToInclude.size() != 3)
        {
            lastError_ = "Error: only three atoms can be specified for angle force calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        int index3 = atoi(atomsToInclude[2].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }
        if((index3 < 0) || (index3 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index3) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFAngleList_.size()))
        {
            lastError_ = std::string("Error: angle force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        C3DVector r3 = state_->atoms_[index3]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();
        std::string atomName3 = state_->atoms_[index3]->getID();

        bool numerical = false;
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "numerical") numerical = true;

        C3DVector f1, f2, f3;
        if(!state_->mdFFAngleList_.get(ffIndex)->isPotentialCalcAvailable() ||
           !state_->mdFFAngleList_.get(ffIndex)->isForceCalcAvailable())
        {
            lastError_ = std::string("Error: angle force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of force calculations!");
            return lastError_;
        }
        if(numerical)
        {
            if(!state_->mdFFAngleList_.get(ffIndex)->calcForcesNumerically(r1, r2, r3, f1, f2, f3))
            {
                printf("Warning: numerical derivatives failed to converge!\r\n");
            }
        }
        else
        {
            state_->mdFFAngleList_.get(ffIndex)->calcForces(r1, r2, r3, f1, f2, f3);
        }

        printf("\r\n");
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f1 = (%g,%g,%g), |f1| = %g\r\n", state_->mdFFAngleList_.get(ffIndex)->getFFType().data(), atomName1.data(), index1, r1.x_, r1.y_, r1.z_, f1.x_, f1.y_, f1.z_, f1.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f2 = (%g,%g,%g), |f2| = %g\r\n", state_->mdFFAngleList_.get(ffIndex)->getFFType().data(), atomName2.data(), index2, r2.x_, r2.y_, r2.z_, f2.x_, f2.y_, f2.z_, f2.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f3 = (%g,%g,%g), |f3| = %g\r\n", state_->mdFFAngleList_.get(ffIndex)->getFFType().data(), atomName3.data(), index3, r3.x_, r3.y_, r3.z_, f3.x_, f3.y_, f3.z_, f3.norm());
        if(numerical) fprintf(stdOut_, "\r\n\tNote: calculations were done numerically!\r\n");
    }
    else if(forceTypeString == "dihedral")
    {
        if(atomsToInclude.size() != 4)
        {
            lastError_ = "Error: only four atoms can be specified for dihedral force calculation!";
            return lastError_;
        }

        int index1 = atoi(atomsToInclude[0].data());
        int index2 = atoi(atomsToInclude[1].data());
        int index3 = atoi(atomsToInclude[2].data());
        int index4 = atoi(atomsToInclude[3].data());
        if((index1 < 0) || (index1 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index1) + std::string(" does not exist!");
            return lastError_;
        }
        if((index2 < 0) || (index2 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index2) + std::string(" does not exist!");
            return lastError_;
        }
        if((index3 < 0) || (index3 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index3) + std::string(" does not exist!");
            return lastError_;
        }
        if((index4 < 0) || (index4 >= state_->atoms_.size()))
        {
            lastError_ = std::string("Error: atom index ") + std::to_string(index4) + std::string(" does not exist!");
            return lastError_;
        }

        text = CASCIIUtility::getArg(arguments, arg++);

        int ffIndex = atoi(text.data());
        if((ffIndex < 0) || (ffIndex >= state_->mdFFDihList_.size()))
        {
            printf("Error: dihedral force field entry %i does not exist!", ffIndex);
            lastError_ = std::string("Error: dihedral force field entry ") + std::to_string(ffIndex) + std::string(" does not exist!");
            return lastError_;
        }

        C3DVector r1 = state_->atoms_[index1]->r_[state_->currentFrame_];
        C3DVector r2 = state_->atoms_[index2]->r_[state_->currentFrame_];
        C3DVector r3 = state_->atoms_[index3]->r_[state_->currentFrame_];
        C3DVector r4 = state_->atoms_[index4]->r_[state_->currentFrame_];
        std::string atomName1 = state_->atoms_[index1]->getID();
        std::string atomName2 = state_->atoms_[index2]->getID();
        std::string atomName3 = state_->atoms_[index3]->getID();
        std::string atomName4 = state_->atoms_[index4]->getID();

        bool numerical = false;
        text = CASCIIUtility::getArg(arguments, arg++);

        if(text == "numerical") numerical = true;

        C3DVector f1, f2, f3, f4;
        if(!state_->mdFFDihList_.get(ffIndex)->isPotentialCalcAvailable() ||
           !state_->mdFFDihList_.get(ffIndex)->isForceCalcAvailable())
        {
            lastError_ = std::string("Error: dihedral force field entry ") + std::to_string(ffIndex) + std::string(" does not have the capability of force calculations!");
            return lastError_;
        }
        if(numerical)
        {
            if(!state_->mdFFDihList_.get(ffIndex)->calcForcesNumerically(r1, r2, r3, r4, f1, f2, f3, f4))
            {
                printf("Warning: numerical derivatives failed to converge!\r\n");
            }
        }
        else
        {
            state_->mdFFDihList_.get(ffIndex)->calcForces(r1, r2, r3, r4, f1, f2, f3, f4);
        }

        printf("\r\n");
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f1 = (%g,%g,%g), |f1| = %g\r\n", state_->mdFFDihList_.get(ffIndex)->getFFType().data(), atomName1.data(), index1, r1.x_, r1.y_, r1.z_, f1.x_, f1.y_, f1.z_, f1.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f2 = (%g,%g,%g), |f2| = %g\r\n", state_->mdFFDihList_.get(ffIndex)->getFFType().data(), atomName2.data(), index2, r2.x_, r2.y_, r2.z_, f2.x_, f2.y_, f2.z_, f2.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f3 = (%g,%g,%g), |f3| = %g\r\n", state_->mdFFDihList_.get(ffIndex)->getFFType().data(), atomName3.data(), index3, r3.x_, r3.y_, r3.z_, f3.x_, f3.y_, f3.z_, f3.norm());
        fprintf(stdOut_, "\tForce (%s) on %s (index %i) at pos (%g,%g,%g): f4 = (%g,%g,%g), |f4| = %g\r\n", state_->mdFFDihList_.get(ffIndex)->getFFType().data(), atomName4.data(), index4, r4.x_, r4.y_, r4.z_, f4.x_, f4.y_, f4.z_, f4.norm());
        if(numerical) fprintf(stdOut_, "\r\n\tNote: calculations were done numerically!\r\n");
    }
    else
    {
        lastError_ = "Error: expected 'nonbonded', 'coulomb', 'bonded', 'angle' or 'dihedral'!";
        return lastError_;
    }

    return lastError_;
}
