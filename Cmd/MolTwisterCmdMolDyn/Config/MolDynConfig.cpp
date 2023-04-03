//
// Copyright (C) 2023 Richard Olsen.
// DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
//
// This file is part of MolTwister.
//
// MolTwister is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MolTwister is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MolTwister.  If not, see <https://www.gnu.org/licenses/>.
//

#include "MolDynConfig.h"

CMolDynConfig::CMolDynConfig()
{
    resetToDefaults();
}

std::vector<std::string> CMolDynConfig::getKeyWords()
{
    return {
             "timestep", "timesteps", "outstride", "ensemble", "NVT", "NPT",
             "temperature", "temperaturerelax", "temperaturenhlen", "temperaturerespa",
             "pressure", "pressurerelax", "pressurenhlen", "pressurerespa",
             "cutoffradius", "neighshell", "cutoffforce",
             "infofile", "xyzfile", "dcdfile", "xtcfile", "pdistrfile", "vdistrfile",
             "includexyzfile", "includedcdfile", "includextcfile", "xtcprecision",
             "includepdistrfile", "includevdistrfile", "maxpdistr", "maxvdistr", "verboseoutput",
             "scale12", "scale13", "scale14", "scale1N", "optmaxsteps", "optlearningrate",
             "optaccuracy"
           };
}

void CMolDynConfig::print(FILE* stdOut)
{
    fprintf(stdOut, "\r\n\tMD configuration parameters\r\n\t------------------------------------------\r\n\r\n");

    fprintf(stdOut, "\ttimesteps = %i; Number of MD timesteps to perform.\r\n", cfg_.numberOfTimeSteps_);
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
    fprintf(stdOut, "\tscale12 = %g; Desired factor to scale 1-2 non-bonded interactions with (i.e., between bonded atoms).\r\n", cfg_.scale12Interactions_);
    fprintf(stdOut, "\tscale13 = %g; Desired factor to scale 1-3 non-bonded interactions with (i.e., between bonded atoms).\r\n", cfg_.scale13Interactions_);
    fprintf(stdOut, "\tscale14 = %g; Desired factor to scale 1-4 non-bonded interactions with (i.e., between bonded atoms).\r\n", cfg_.scale14Interactions_);
    fprintf(stdOut, "\tscale1N = %g; Desired factor to scale 1-N>4 non-bonded interactions with (i.e., between bonded atoms).\r\n", cfg_.scaleAbove14BondedInteractions_);
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\tinfofile = %s; Filename of info file.\r\n", cfg_.outInfoFile_.data());
    fprintf(stdOut, "\txyzfile = %s; Filename of XYZ file.\r\n", cfg_.outXYZFile_.data());
    fprintf(stdOut, "\tdcdfile = %s; Filename of DCD file.\r\n", cfg_.outDCDFile_.data());
    fprintf(stdOut, "\txtcfile = %s; Filename of Gromacs XTC file.\r\n", cfg_.outXTCFile_.data());
    fprintf(stdOut, "\tpdistrfile = %s; Filename of momentum distribution file.\r\n", cfg_.outPDistrFile_.data());
    fprintf(stdOut, "\tvdistrfile = %s; Filename of volume distribution file.\r\n", cfg_.outVDistrFile_.data());
    fprintf(stdOut, "\tincludexyzfile {yes, no} = %s; Include XYZ file as output.\r\n", cfg_.includeXYZFile_ ? "Yes" : "No");
    fprintf(stdOut, "\tincludedcdfile {yes, no} = %s; Include DCD file as output.\r\n", cfg_.includeDCDFile_ ? "Yes" : "No");
    fprintf(stdOut, "\tincludextcfile {yes, no} = %s; Include Gromacs XTC file as output.\r\n", cfg_.includeXTCFile_ ? "Yes" : "No");
    fprintf(stdOut, "\txtcprecision = %g; Gromacs XTC file precision.\r\n", cfg_.xtcPrecision_);
    fprintf(stdOut, "\tincludepdistrfile {yes, no} = %s; Include momentum distribution file as output.\r\n", cfg_.includePDistrFile_ ? "Yes" : "No");
    fprintf(stdOut, "\tincludevdistrfile {yes, no} = %s; Include volume distribution file as output.\r\n", cfg_.includeVDistrFile_ ? "Yes" : "No");
    fprintf(stdOut, "\tmaxpdistr = %g AA*g/(fs*mol); Desired maximum momentum to use for momentum distribution output.\r\n", cfg_.maxPDistrOutput_);
    fprintf(stdOut, "\tmaxvdistr = %g AA^3; Desired maximum volume to use for volume distribution output.\r\n", cfg_.maxVDistrOutput_);
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\tverboseoutput {yes, no} = %s; Let output to screen be verbose.\r\n", cfg_.verboseOutput_ ? "Yes" : "No");
    fprintf(stdOut, "\r\n");
    fprintf(stdOut, "\toptlearningrate = %g; Learning rate, or step size, of gradient descent algorithm used for energy optimization.\r\n", cfg_.gradientDescentLearningRate_);
    fprintf(stdOut, "\toptmaxsteps = %i; Maximum number of steps used for energy optimization.\r\n", cfg_.gradientDescentMaxSteps_);
    fprintf(stdOut, "\toptaccuracy = %g; Accuracy, a, for energy optimization in kJ/mol, which terminates if |U_(n+1) - U_(n)| < a.\r\n", cfg_.gradientDescentAccuracy_);

    fprintf(stdOut, "\r\n\t------------------------------------------\r\n");
}

std::string CMolDynConfig::set(std::string parameter, std::string value)
{
    if(parameter == "timesteps")
    {
        cfg_.numberOfTimeSteps_ = std::atoi(value.data());
    }
    else if(parameter == "timestep")
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
    else if(parameter == "dcdfile")
    {
        cfg_.outDCDFile_ = value.data();
    }
    else if(parameter == "xtcfile")
    {
        cfg_.outXTCFile_ = value.data();
    }
    else if(parameter == "pdistrfile")
    {
        cfg_.outPDistrFile_ = value.data();
    }
    else if(parameter == "vdistrfile")
    {
        cfg_.outVDistrFile_ = value.data();
    }
    else if(parameter == "includexyzfile")
    {
        if(value == "yes") cfg_.includeXYZFile_ = true;
        else if(value == "no") cfg_.includeXYZFile_ = false;
        else
        {
            return "Error: 'includexyzfile' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "includedcdfile")
    {
        if(value == "yes") cfg_.includeDCDFile_ = true;
        else if(value == "no") cfg_.includeDCDFile_ = false;
        else
        {
            return "Error: 'includedcdfile' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "includextcfile")
    {
        if(value == "yes") cfg_.includeXTCFile_ = true;
        else if(value == "no") cfg_.includeXTCFile_ = false;
        else
        {
            return "Error: 'includextcfile' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "xtcprecision")
    {
        cfg_.xtcPrecision_ = std::atof(value.data());
    }
    else if(parameter == "includepdistrfile")
    {
        if(value == "yes") cfg_.includePDistrFile_ = true;
        else if(value == "no") cfg_.includePDistrFile_ = false;
        else
        {
            return "Error: 'includepdistrfile' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "includevdistrfile")
    {
        if(value == "yes") cfg_.includeVDistrFile_ = true;
        else if(value == "no") cfg_.includeVDistrFile_ = false;
        else
        {
            return "Error: 'includevdistrfile' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "maxpdistr")
    {
        cfg_.maxPDistrOutput_ = std::atof(value.data());
    }
    else if(parameter == "maxvdistr")
    {
        cfg_.maxVDistrOutput_ = std::atof(value.data());
    }
    else if(parameter == "verboseoutput")
    {
        if(value == "yes") cfg_.verboseOutput_ = true;
        else if(value == "no") cfg_.verboseOutput_ = false;
        else
        {
            return "Error: 'verboseoutput' parameter should be either 'yes' or 'no'!";
        }
    }
    else if(parameter == "scale12")
    {
        cfg_.scale12Interactions_ = std::atof(value.data());
    }
    else if(parameter == "scale13")
    {
        cfg_.scale13Interactions_ = std::atof(value.data());
    }
    else if(parameter == "scale14")
    {
        cfg_.scale14Interactions_ = std::atof(value.data());
    }
    else if(parameter == "scale1N")
    {
        cfg_.scaleAbove14BondedInteractions_ = std::atof(value.data());
    }
    else if(parameter == "optlearningrate")
    {
        cfg_.gradientDescentLearningRate_ = std::atof(value.data());
    }
    else if(parameter == "optmaxsteps")
    {
        cfg_.gradientDescentMaxSteps_ = std::atoi(value.data());
    }
    else if(parameter == "optaccuracy")
    {
        cfg_.gradientDescentAccuracy_ = std::atof(value.data());
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
    cfg_.cutoffRadius_ = 10.0; // AA
    cfg_.neighListShell_ = 2.0; // AA
    cfg_.cutoffForce_ = 1000.0; // Reduced units
    cfg_.numberOfTimeSteps_ = 50000;
    cfg_.outInfoFile_ = "out.txt";
    cfg_.outXYZFile_ = "traj.xyz";
    cfg_.outDCDFile_ = "traj.dcd";
    cfg_.outXTCFile_ = "traj.xtc";
    cfg_.outPDistrFile_ = "distr";
    cfg_.outVDistrFile_ = "distr";
    cfg_.includeXYZFile_ = false;
    cfg_.includeDCDFile_ = true;
    cfg_.includeXTCFile_ = false;
    cfg_.xtcPrecision_ = 1000.0f;
    cfg_.includePDistrFile_ = false;
    cfg_.includeVDistrFile_ = false;
    cfg_.maxPDistrOutput_ = 0.3; // AA*g/(fs*mol)
    cfg_.maxVDistrOutput_ = 100.0E3; // AA^3
    cfg_.verboseOutput_ = false;
    cfg_.scale12Interactions_ = 0.0;
    cfg_.scale13Interactions_ = 0.0;
    cfg_.scale14Interactions_ = 0.5;
    cfg_.scaleAbove14BondedInteractions_ = 0.0;
    cfg_.gradientDescentLearningRate_ = 1e-6;
    cfg_.gradientDescentMaxSteps_ = 5000;
    cfg_.gradientDescentAccuracy_ = 0.0;
}
