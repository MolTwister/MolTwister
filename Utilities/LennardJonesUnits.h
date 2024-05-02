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

class CLJUnits
{
public:
    CLJUnits() { defineUnits(); }
    
public:
    double massUnit() const { return massUnit_; }
    double distanceUnit() const { return distanceUnit_; }
    double energyUnit() const { return energyUnit_; }
    double timeUnit() const { return timeUnit_; }
    double velocityUnit() const { return velocityUnit_; }
    double forceUnit() const { return forceUnit_; }
    double torqueUnit() const { return torqueUnit_; }
    double tempUnit() const { return tempUnit_; }
    double pressUnit() const { return pressUnit_; }
    double chargeUnit() const { return chargeUnit_; }
    double volumeUnit() const { return volumeUnit_; }
    double buckinghamC() const { return buckinghamC_; }
    double harmonicBondK() const { return harmonicBondK_; }
    
private:
    double massUnit_;
    double distanceUnit_;
    double energyUnit_;
    double timeUnit_;
    double velocityUnit_;
    double forceUnit_;
    double torqueUnit_;
    double tempUnit_;
    double pressUnit_;
    double chargeUnit_;
    double volumeUnit_;
    double buckinghamC_;
    double harmonicBondK_;

private:
    void defineUnits();
};
