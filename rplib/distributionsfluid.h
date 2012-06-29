#ifndef RPLIB_DISTRIBUTIONS_FLUID_H
#define RPLIB_DISTRIBUTIONS_FLUID_H

#include "rplib/fluid.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluid {
public:

  DistributionsFluid(){}

  virtual ~DistributionsFluid(){}

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid * GenerateSample(const std::vector<double> & /*trend_params*/) const = 0;

};

#endif
