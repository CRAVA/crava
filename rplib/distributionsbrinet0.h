#ifndef DISTRIBUTIONSBRINET0_H
#define DISTRIBUTIONSBRINET0_H

#include "rplib/fluid.h"

// Parallel classes are Brine and DistributionsBrineEvolution.
class DistributionsBrineT0 : public DistributionsFluidT0 {
public:

  DistributionsBrineT0() : DistributionsFluidT0(){}

  virtual ~DistributionsBrineT0(){}

  virtual Fluid * GenerateSample() const {
    double temperature   = 30.0;     //FAKE
    double pore_pressure = 10.0;     //FAKE
    double salinity      = 0.05;     //FAKE
    Fluid * fluid = new Brine(temperature, pore_pressure, salinity);
    return fluid;
  }
};

#endif
