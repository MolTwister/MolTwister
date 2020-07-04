#pragma once
#include "MolTwisterCmd.h"

class CCmdMolDyn : public CCmd
{
public:
    class CConfig
    {
    public:
        enum Ensemble { ensembleNVT=0, ensembleNPT=1 };

    public:
        CConfig()
        {
            timeStep_ = 1.0; // fs
            ensemble_ = ensembleNVT;
            temperature_ = 298.0; // K
            temperatureRelaxationTime_ = 20.0 * timeStep_; // fs
            temperatureNHChainLength_ = 4;
            temperatureRESPASteps_ = 4;
            pressure_ = 1.0; // atm
            pressureRelaxationTime_ = 5000.0 * timeStep_; // fs
            pressureNHChainLength_ = 4;
            pressureRESPASteps_ = 4;
            cutoffRadius_ = 10.0; // Å
            cuoffForce_ = 1000.0; // Reduced units
            numberOfTimeSteps_ = 10000;
        }

    public:
        void print(FILE* stdOut);

    public:
        Ensemble ensemble_;
        double temperature_;
        double temperatureRelaxationTime_;
        int temperatureNHChainLength_;
        int temperatureRESPASteps_;
        double pressure_;
        double pressureRelaxationTime_;
        int pressureNHChainLength_;
        int pressureRESPASteps_;
        double cutoffRadius_;
        double cuoffForce_;
        int numberOfTimeSteps_;
        double timeStep_;
    };

public:
    CCmdMolDyn() = delete;
    CCmdMolDyn(CMolTwisterState* state) : CCmd(state) { state_ = state; }

public:
    void execute(std::string commandLine);
    std::string getCmd() const { return "moldyn"; }
    std::string getTopLevHelpString() const { return std::string("Invoke a molecular dynamics simulation"); }
    std::string getHelpString() const;
    std::string getHelpString(std::string subCommand) const;

protected:
    virtual void onAddKeywords();
    virtual void onRegisterSubCommands();

private:
    CConfig config_;
};
