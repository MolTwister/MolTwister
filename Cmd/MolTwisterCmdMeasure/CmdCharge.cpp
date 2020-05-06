#include "CmdCharge.h"
#include "../../Utilities/ASCIIUtility.h"

namespace MolTwisterCmdMeasure
{
    std::string CCmdCharge::getCmd()
    {
        return "charge";
    }

    std::vector<std::string> CCmdCharge::getCmdLineKeywords()
    {
        return { "charge", "sel", "tot" };
    }

    std::vector<std::string> CCmdCharge::getCmdHelpLines()
    {
        return {
                    "charge sel",
                    "charge tot"
               };
    }

    std::string CCmdCharge::getCmdFreetextHelp()
    {
        std::string text;

        text+= "\tMeasures the total charge of the collection of atoms defined by\r\n";
        text+= "\t* the visual selection of atoms, by using the 'sel' keyword\r\n";
        text+= "\t* all the atoms, by using the 'tot' keyword\r\n";
        text+= "\tIf 'sel' or 'tot' is not specified, 'tot' is assumed.\r\n";
        text+= "\r\n";
        text+= "\tOutput:\r\n";
        text+= "\tQsum = <summed charge>";

        return text;
    }

    std::string CCmdCharge::execute(std::vector<std::string> arguments)
    {
        lastError_ = "";

        size_t arg = 0;
        std::string text;
        double Q = 0.0;
        int typeOfSelection = 0;

        // Type of selections possible are
        // 0: of total (default), 1: of selected atoms
        text = CASCIIUtility::getArg(arguments, arg++);
        if(text.size() == 0)    typeOfSelection = 0;
        else if(text == "sel")  typeOfSelection = 1;
        else if(text == "tot")  typeOfSelection = 0;
        else
        {
            lastError_ = "Error, type of selection not recognized!";
            return lastError_;
        }

        if(typeOfSelection == 1)
        {
            for(int i=0; i<state_->atoms_.size(); i++)
            {
                if(state_->atoms_[i]->isSelected())
                    Q+= state_->atoms_[i]->Q_;
            }
        }
        else
        {
            for(int i=0; i<state_->atoms_.size(); i++)
            {
                Q+= state_->atoms_[i]->Q_;
            }
        }

        fprintf(stdOut_, "\r\n\tQsum = %.8f\r\n", Q);

        return lastError_;
    }
}
