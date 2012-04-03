#ifndef RPLIB_DISTRIBUTIONSCO2EVOLUTION_H
#define RPLIB_DISTRIBUTIONSCO2EVOLUTION_H

#include "rplib/distributionsfluidevolution.h"

// Parallel classes are CO2 and DistributionsCO2.
class DistributionsCO2Evolution : public DistributionsFluidEvolution {
public:

  DistributionsCO2Evolution() : DistributionsFluidEvolution(){}

  virtual ~DistributionsCO2Evolution(){}

  double Sample() const {return 1.0;}  // FAKE

};

#endif
