#ifndef DISTRIBUTIONSBRINE_H
#define DISTRIBUTIONSBRINE_H

#include "rplib/distributionsfluid.h"
#include "rplib/brine.h"

// Parallel classes are Brine and DistributionsBrineEvolution.
class DistributionsBrine : public DistributionsFluid {
public:

  DistributionsBrine() : DistributionsFluid(){}

  virtual ~DistributionsBrine(){}

  virtual Fluid * GenerateSample() const {
    double salinity      = 0.05;     //FAKE
    Fluid * fluid = new Brine(salinity);
    return fluid;
  }
};

#endif
