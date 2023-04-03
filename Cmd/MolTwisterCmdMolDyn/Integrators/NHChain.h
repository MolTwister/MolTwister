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

#pragma once
#include <stdio.h>
#include <vector>
#include "../../../Utilities/CUDAGeneralizations.h"

BEGIN_CUDA_COMPATIBLE()

class CFct
{
public:
    CFct() {}
    virtual ~CFct() {}
    
public:
    virtual double G(int j, std::vector<double>& p_eta, std::vector<double>& Q, double beta) = 0;
    virtual void scaleMomentum(double coeff) = 0;
};

class CNHChain
{
public:
    CNHChain();
    
public:
    void propagator(int N, int dim, double dt, CFct& f);
    void prepareArrays(int N, int dim);
    void setRandNHPos();
    void setRandNHMom();
    void setVerboseOutput(bool verbose=true) { verboseOutput_ = verbose; }

private:
    void getSuzukiYoshida(double* w_);
    
public:
    double T_;                                       // Temperature in reduced units
    int n_;                                          // RESPA steps in Nose-Hoover part
    int M_;                                          // Nose-Hoover chain length
    double tau_;                                     // Thermostat relaxation time in reduced units
    std::vector<double> eta_;                        // Nose-Hoover position coordinates
    std::vector<double> p_eta_;                      // Nose-Hoover momentum
    std::vector<double> Q_;                          // Nose-Hoover Q for each chain entry
    
private:
    int n_sy_;
    double w_[3];
    bool verboseOutput_;
};

END_CUDA_COMPATIBLE()
