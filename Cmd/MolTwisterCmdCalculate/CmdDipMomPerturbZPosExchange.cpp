#include "CmdDipMomPerturbZPosExchange.h"
#include "../../Utilities/ASCIIUtility.h"
#include "../Tools/DCDTools.h"
#include "../Tools/MolTwisterStateTools.h"

std::string CCmdDipMomPerturbZPosExchange::getCmd()
{
    return "dipmomperturbzposexchange";
}

std::vector<std::string> CCmdDipMomPerturbZPosExchange::getCmdLineKeywords()
{
    return { "dipmomperturbzposexchange", "chargedmol", "__fromloadedframes__" };
}

std::vector<std::string> CCmdDipMomPerturbZPosExchange::getCmdHelpLines()
{
    return {
                "dipmomperturbzposexchange <DCD filename> <frame index> <list of atomic IDs> <positive charges> <negative charges> <delta Z pos. charges> <delta Z neg. charges> [chargedmol]"
           };
}

std::string CCmdDipMomPerturbZPosExchange::getCmdFreetextHelp()
{
    std::string text;

    text+= "\tCalculate dipole moment while incrementally increasing the number of ion-swaps to\r\n";
    text+= "\tperform. Swapping the ions will naturally keep the system charge neutral. The ion\r\n";
    text+= "\tindices to exchange are chosen at random. I.e., in the first iteration a random\r\n";
    text+= "\tion-pair is chosen. Then, in the next iteration, the first pair remains the same\r\n";
    text+= "\twhile the second pair is randomly chosen. An ion pair is constructed by first\r\n";
    text+= "\tselecting an index from the positive ions and then finding the negative ion\r\n";
    text+= "\twhich is the closest to the chosen positive ion. The the dipole moment is calculated\r\n";
    text+= "\tfor the selected molecules <list of atomic IDs> (comma separated, no space). The\r\n";
    text+= "\tdipolemoment is averaged based on all the defined molecules in frame index given by\r\n";
    text+= "\t<frame index>. The DCD file, <DCD filename>, is used as input, or the loaded frames if\r\n";
    text+= "\t<DCD filename> = __fromloadedframes__. The dipole moment expression for neutral molecules\r\n";
    text+= "\tare used as default. By sepcifying 'chargedmol', the dipole moment expression for charged\r\n";
    text+= "\tmolecules is employed. The lists of positive and negative charges are lists of atomic IDs,\r\n";
    text+= "\tsuch as H, O and C7 (comma separated, no space), which defines groups of positive and\r\n";
    text+= "\tnegative charges that can form the aforementioned ion pairs to swap (i.e., ions are picked\r\n";
    text+= "\tfrom these groups). The <delta Z pos. charges> and <delta Z neg. charges> parameters\r\n";
    text+= "\tdefine the z-axis displacement to perform on the positive and negative charges,\r\n";
    text+= "\trespectively.\r\n";
    text+= "\r\n";
    text+= "\tOutput:\r\n";
    text+= "\t1. #Exchanges Px Py Pz Exch.index\r\n";
    text+= "\t2. <number of performed exchanges> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Index of the positive ion>\r\n";
    text+= "\t3. <number of performed exchanges> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Index of the positive ion>\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\t                 .\r\n";
    text+= "\tN+1. <number of performed exchanges> <dipole moment x-component> <dipole moment y-component> <dipole moment z-component> <Index of the positive ion>\r\n";
    text+= "\tN+2.             .\r\n";
    text+= "\tN+3.             .\r\n";
    text+= "\t             <N atomic system configurations, for each exchange step that has been performed, in the xyz-file format>\r\n";
    text+= "\tM.               .\r\n";
    text+= "\tM+1.             .\r\n";
    text+= "\twhere N is the number of perturbations.";

    return text;
}

std::string CCmdDipMomPerturbZPosExchange::execute(std::vector<std::string> arguments)
{
    lastError_ = "";

    size_t arg = 0;
    std::string text;
    CDCDFile dcdFile;

    // Open DCD file
    text = CASCIIUtility::getArg(arguments, arg++);

    bool useDCD = false;
    if(text != "__fromloadedframes__")
    {
        if(!dcdFile.open(text))
        {
            lastError_ = std::string("Error: Unable to open file ") + text + std::string("!");
            return lastError_;
        }

        useDCD = true;
    }

    text = CASCIIUtility::getArg(arguments, arg++);
    int frame = atoi(text.data());

    if(useDCD)
    {
        if((frame < 0) || (frame >= dcdFile.getNumRecords()))
        {
            lastError_ = "Error: Frame is outside the range of available number of frames (note: index should be zero indexed)";
            return lastError_;
        }
    }
    else
    {
        if(state_->atoms_.size() < 1)
        {
            lastError_ = "Error: no atoms found";
            return lastError_;
        }
        if((frame < 0) || (frame >= state_->atoms_[0]->r_.size()))
        {
            lastError_ = "Error: Frame is outside the range of available number of frames (note: index should be zero indexed)";
            return lastError_;
        }
    }

    // Get atom types to calculate the dipole moment for
    std::vector<std::string> atomsDipmomBase;
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsDipmomBase = CASCIIUtility::getWords(text, ",");

    // Get atom types with positive charge
    std::vector<std::string> atomsPositive;
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsPositive = CASCIIUtility::getWords(text, ",");

    // Get atom types with negative charge
    std::vector<std::string> atomsNegative;
    text = CASCIIUtility::getArg(arguments, arg++);
    CASCIIUtility::removeWhiteSpace(text);
    atomsNegative = CASCIIUtility::getWords(text, ",");

    // Get the z-distance to displace the positive and negative charges, respectively
    text = CASCIIUtility::getArg(arguments, arg++);
    double displacementPositives = atof(text.data());

    text = CASCIIUtility::getArg(arguments, arg++);
    double displacementNegatives = atof(text.data());

    // Are we going to use charge neutral formulation: sum_i q_i r_i or
    // are we using charged molecules: sum_i q_i (r_i - R_c)?
    bool chargeNeutralFormulation = true;
    text = CASCIIUtility::getArg(arguments, arg++);

    if(text == "chargedmol")
        chargeNeutralFormulation = false;
    else arg--;

    // List all selected atom indices to calculate the dipole moment for
    std::vector<int> indicesDipmomBase;
    for(int i=0; i<(int)atomsDipmomBase.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsDipmomBase[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesDipmomBase.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // List all selected atom indices associated with a positive charge
    std::vector<int> indicesPositive;
    for(int i=0; i<(int)atomsPositive.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsPositive[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesPositive.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // List all selected atom indices associated with a negative charge
    std::vector<int> indicesNegative;
    for(int i=0; i<(int)atomsNegative.size(); i++)
    {
        std::vector<CAtom*> atoms;
        state_->getAtomsWithID(atomsNegative[i], atoms);
        for(int j=0; j<(int)atoms.size(); j++)
        {
            indicesNegative.emplace_back(atoms[j]->getAtomIndex());
        }
    }

    // Retrieve DCD at selected frame
    if(useDCD) dcdFile.gotoRecord(frame);

    // Print header
    fprintf(stdOut_, "%-20s%-20s%-20s%-20s%-20s\r\n", "#Exchanges", "Px", "Py", "Pz", "Exch.index");

    // Print the dipolemoment without any perturbations
    C3DVector P, Rc;
    CDCDTools dcdTools(state_, stdOut_);
    if(useDCD) P = dcdTools.getMoleculeDipoleMoment(indicesDipmomBase, &dcdFile, Rc, chargeNeutralFormulation);
    else P = CMolTwisterStateTools(state_, stdOut_).getMoleculeDipoleMoment(indicesDipmomBase, frame, Rc, chargeNeutralFormulation);
    fprintf(stdOut_, "%-20i%-20.8f%-20.8f%-20.8f%-20i\r\n", 0, P.x_, P.y_, P.z_, -1);

    // Calculate dipole moment while incrementally increasing the number of ion-swaps to
    // perform. Swapping the ions will naturally keep the system charge neutral. The ion
    // indices to exchange are chosen at random. I.e., in the first iteration a random
    // ion-pair is chosen. Then, in the next iteration, the first pair remains the same
    // while the second pair is randomly chosen. An ion pair is constructed by first
    // selecting an index from the positive ions and then finding the negative ion
    // which is the closest to the chosen positive ion.
    srand(static_cast<unsigned int>(time(nullptr)));
    int lowestIonCount = static_cast<int>((indicesPositive.size() < indicesNegative.size()) ? indicesPositive.size() : indicesNegative.size());
    std::vector<int> selectedIonIndices;
    std::vector<std::vector<std::pair<std::string, C3DVector> > > configurations(lowestIonCount);
    for(int i=0; i<lowestIonCount; i++)
    {
        // Select new index until finding one that does not exist
        int ionIndex = -1;
        bool doesExist = false;
        do
        {
            ionIndex = static_cast<int>((static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) * static_cast<double>(lowestIonCount));
            doesExist = false;
            for(int j=0; j<selectedIonIndices.size(); j++)
            {
                if(ionIndex == selectedIonIndices[j])
                {
                    doesExist = true;
                    break;
                }
            }

            if(!doesExist) selectedIonIndices.emplace_back(ionIndex);

        } while(doesExist);

        // Get position of positive ion
        C3DVector r_p;
        if(useDCD) r_p = dcdFile.getCoordinate(indicesPositive[ionIndex]);
        else r_p = state_->atoms_[indicesPositive[ionIndex]]->r_[frame];

        // Find the negative ion, which is closest (in the x/y-plane) to the chosen positive ion
        int associatedNegativeIon = -1;
        double closestDistance2 = DBL_MAX;
        for(int j=0; j<indicesNegative.size(); j++)
        {
            C3DVector r_n;
            if(useDCD) r_n = dcdFile.getCoordinate(indicesNegative[j]);
            else r_n = state_->atoms_[indicesNegative[j]]->r_[frame];
            double dx = r_n.x_ - r_p.x_;
            double dy = r_n.y_ - r_p.y_;
            double dist = dx*dx + dy*dy;
            if(dist < closestDistance2)
            {
                closestDistance2 = dist;
                associatedNegativeIon = j;
            }
        }

        // We have found an ion pair. Now we need to exchange the ions.
        C3DVector r_n;
        if(useDCD) r_n = dcdFile.getCoordinate(indicesNegative[associatedNegativeIon]);
        else r_n = state_->atoms_[indicesNegative[associatedNegativeIon]]->r_[frame];

        r_p.z_+= displacementPositives;
        r_n.z_+= displacementNegatives;
        if(useDCD)
        {
            dcdFile.setCoordinate(indicesPositive[ionIndex], r_p);
            dcdFile.setCoordinate(indicesNegative[associatedNegativeIon], r_n);
        }
        else
        {
            state_->atoms_[indicesPositive[ionIndex]]->r_[frame] = r_p;
            state_->atoms_[indicesNegative[associatedNegativeIon]]->r_[frame] = r_n;
        }

        // Measure the dipole moment and print it
        C3DVector P, Rc;
        if(useDCD) P = dcdTools.getMoleculeDipoleMoment(indicesDipmomBase, &dcdFile, Rc, chargeNeutralFormulation);
        else P = CMolTwisterStateTools(state_, stdOut_).getMoleculeDipoleMoment(indicesDipmomBase, frame, Rc, chargeNeutralFormulation);
        fprintf(stdOut_, "%-20i%-20.8f%-20.8f%-20.8f%-20i\r\n", i+1, P.x_, P.y_, P.z_, indicesPositive[ionIndex]);

        // Store the current configuration for later printing
        int numCoordinates;
        if(useDCD) numCoordinates = dcdFile.getNumCoordinatesInRecord();
        else numCoordinates = static_cast<int>(state_->atoms_.size());
        std::vector<std::pair<std::string, C3DVector> > configuration(numCoordinates);
        for(int j=0; j<numCoordinates; j++)
        {
            configuration[j].first = state_->atoms_[j]->getID();
            if(useDCD) configuration[j].second = dcdFile.getCoordinate(j);
            else configuration[j].second = state_->atoms_[j]->r_[frame];
        }
        configurations[i] = configuration;
    }

    // Print configurations for each step
    fprintf(stdOut_, "\r\n");
    for(int i=0; i<configurations.size(); i++)
    {
        fprintf(stdOut_, "%i\r\n", static_cast<int>(configurations[i].size()));
        fprintf(stdOut_, "Number of exchanges made = %i\r\n", i+1);
        for(int j=0; j<configurations[i].size(); j++)
        {
            fprintf(stdOut_, "%-20s% -20.8f% -20.8f% -20.8f\r\n", configurations[i][j].first.data(), configurations[i][j].second.x_, configurations[i][j].second.y_, configurations[i][j].second.z_);
        }
        fprintf(stdOut_, "\r\n");
    }

    return lastError_;
}
