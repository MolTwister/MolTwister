#include "MolDynConfig.h"

CMolDynConfig::CMolDynConfig()
{
    resetToDefaults();
}

std::vector<std::string> CMolDynConfig::getKeyWords()
{
    return {
             "timestep", "outstride", "ensemble", "NVT", "NPT",
             "temperature", "temperaturerelax", "temperaturenhlen", "temperaturerespa",
             "pressure", "pressurerelax", "pressurenhlen", "pressurerespa",
        "cutoffradius", "neighshell", "cutoffforce",
             "infofile", "xyzfile"
           };
}

void CMolDynConfig::print(FILE* stdOut)
{
    fprintf(stdOut, "\r\n\tMD configuration parameters\r\n\t------------------------------------------\r\n\r\n");

    fprintf(stdOut, "\ttimestep = %g fs; Time step between each MD integration.\r\n", cfg_.timeStep_);
    fprintf(stdOut, "\toutstride = %i; Number of time steps between each trajectory output.\r\n", cfg_.outputStride_);
    fprintf(stdOut, "\tensemble {nve, nvt, npt} = %s; Ensemble, where n, v, p, e and t are atom count, volume, pressure, energy and temperature, respectivley.\r\n", (cfg_.ensemble_ == SMolDynConfigStruct::ensembleNVT) ? "NVT" : ((cfg_.ensemble_ == SMolDynConfigStruct::ensembleNPT) ? "NPT" : (cfg_.ensemble_ == SMolDynConfigStruct::ensembleNVE) ? "NVE" : "Unknown"));
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\ttemperature = %g K; Desired temperature of system.\r\n", cfg_.temperature_);
    fprintf(stdOut, "\ttemperaturerelax = %g fs; Nose-Hoover thermal relaxation time.\r\n", cfg_.temperatureRelaxationTime_);
    fprintf(stdOut, "\ttemperaturenhlen = %i; Nose-Hoover thermal chain length.\r\n", cfg_.temperatureNHChainLength_);
    fprintf(stdOut, "\ttemperaturerespa = %i; Nose-Hoover thermal RESPA steps.\r\n", cfg_.temperatureRESPASteps_);
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\tpressure = %g atm; Desired pressure of system.\r\n", cfg_.pressure_);
    fprintf(stdOut, "\tpressurerelax = %g fs; Nose-Hoover pressure relaxation time.\r\n", cfg_.pressureRelaxationTime_);
    fprintf(stdOut, "\tpressurenhlen = %i; Nose-Hoover pressure chain length.\r\n", cfg_.pressureNHChainLength_);
    fprintf(stdOut, "\tpressurerespa = %i; Nose-Hoover pressure RESPA steps.\r\n", cfg_.pressureRESPASteps_);
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\tcutoffradius = %g AA; Desired cutoff radius.\r\n", cfg_.cutoffRadius_);
    fprintf(stdOut, "\tneighshell = %g AA; Desired neighbor list shell distance.\r\n", cfg_.neighListShell_);
    fprintf(stdOut, "\tcutoffforce = %g; Desired cutoff force (in reduced units).\r\n", cfg_.cutoffForce_);
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\tinfofile = %s; Filename of info file.\r\n", cfg_.outInfoFile_.data());
    fprintf(stdOut, "\txyzfile = %s; Filename of XYZ file.\r\n", cfg_.outXYZFile_.data());

    fprintf(stdOut, "\r\n\t------------------------------------------\r\n");
}

std::string CMolDynConfig::set(std::string parameter, std::string value)
{
    if(parameter == "timestep")
    {
        cfg_.timeStep_ = std::atof(value.data());
    }
    else if(parameter == "outstride")
    {
        cfg_.outputStride_ = std::atoi(value.data());
    }
    else if(parameter == "ensemble")
    {
        if(value == "nvt") cfg_.ensemble_ = SMolDynConfigStruct::ensembleNVT;
        else if(value == "npt") cfg_.ensemble_ = SMolDynConfigStruct::ensembleNPT;
        else if(value == "nve") cfg_.ensemble_ = SMolDynConfigStruct::ensembleNVE;
        else
        {
            return "Error: unknown ensemble!";
        }
    }
    else if(parameter == "temperature")
    {
        cfg_.temperature_ = std::atof(value.data());
    }
    else if(parameter == "temperaturerelax")
    {
        cfg_.temperatureRelaxationTime_ = std::atof(value.data());
    }
    else if(parameter == "temperaturenhlen")
    {
        cfg_.temperatureNHChainLength_ = std::atoi(value.data());
    }
    else if(parameter == "temperaturerespa")
    {
        cfg_.temperatureRESPASteps_ = std::atoi(value.data());
    }
    else if(parameter == "pressure")
    {
        cfg_.pressure_ = std::atof(value.data());
    }
    else if(parameter == "pressurerelax")
    {
        cfg_.pressureRelaxationTime_ = std::atof(value.data());
    }
    else if(parameter == "pressurenhlen")
    {
        cfg_.pressureNHChainLength_ = std::atoi(value.data());
    }
    else if(parameter == "pressurerespa")
    {
        cfg_.pressureRESPASteps_ = std::atoi(value.data());
    }
    else if(parameter == "cutoffradius")
    {
        cfg_.cutoffRadius_ = std::atof(value.data());
    }
    else if(parameter == "neighshell")
    {
        cfg_.neighListShell_ = std::atof(value.data());
    }
    else if(parameter == "cutoffforce")
    {
        cfg_.cutoffForce_ = std::atof(value.data());
    }
    else if(parameter == "infofile")
    {
        cfg_.outInfoFile_ = value.data();
    }
    else if(parameter == "xyzfile")
    {
        cfg_.outXYZFile_ = value.data();
    }
    else
    {
        return "Error: unknown parameter!";
    }

    return "";
}

void CMolDynConfig::resetToDefaults()
{
    cfg_.timeStep_ = 1.0; // fs
    cfg_.outputStride_ = 100;
    cfg_.ensemble_ = cfg_.ensembleNVT;
    cfg_.temperature_ = 298.0; // K
    cfg_.temperatureRelaxationTime_ = 20.0 * cfg_.timeStep_; // fs
    cfg_.temperatureNHChainLength_ = 4;
    cfg_.temperatureRESPASteps_ = 4;
    cfg_.pressure_ = 1.0; // atm
    cfg_.pressureRelaxationTime_ = 5000.0 * cfg_.timeStep_; // fs
    cfg_.pressureNHChainLength_ = 4;
    cfg_.pressureRESPASteps_ = 4;
    cfg_.cutoffRadius_ = 10.0; // Å
    cfg_.neighListShell_ = 2.0; // Å
    cfg_.cutoffForce_ = 1000.0; // Reduced units
    cfg_.numberOfTimeSteps_ = 50000;
    cfg_.outInfoFile_ = "out.txt";
    cfg_.outXYZFile_ = "traj.xyz";
}
