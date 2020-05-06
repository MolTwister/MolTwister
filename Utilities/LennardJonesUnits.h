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
