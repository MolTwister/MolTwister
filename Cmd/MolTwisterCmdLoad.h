#pragma once
#include "MolTwisterCmd.h"

class CCmdLoad : public CCmd
{
public:
    CCmdLoad() = delete;
    CCmdLoad(CMolTwisterState* state) : CCmd(state) { state_ = state; }
    
public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "load"; }
    std::string getTopLevHelpString() const { return std::string("Load data into MolTwister"); }
    std::string getHelpString() const;
    
protected:
    virtual void onAddKeywords();

private:
    void parseXYZCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parsePDBCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseMTTCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseScriptCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parsePythonCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseMasschargeCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);
    void parseQeposCommand(std::string commandLine, int& arg, std::vector<std::string>& bondAtomsToIgnore, bool& genBonds, bool& ignoreAllBonds, int& baseFrameIndex);

    bool readXYZFile(std::string xyzFileName, bool& genBonds);
    bool readPDBFile(std::string pdbFileName, bool& genBonds);
    bool readMTTFile(std::string mttFileName);
    bool readScriptFile(std::string scriptFileName);
    bool readPythonFile(std::string scriptFileName);
    bool readMassChargeFile(std::string massChargeFileName);
    bool readQEPosFile(std::string qePosFileName, std::string qeInputFileName);
};
