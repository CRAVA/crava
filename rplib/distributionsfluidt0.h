#ifndef DISTRIBUTIONSFLUIDT0_H
#define DISTRIBUTIONSFLUIDT0_H

#include "rplib/fluid.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluidT0 {
public:

  DistributionsFluidT0(){}

  virtual ~DistributionsFluidT0(){}

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid * GenerateSample() const = 0;
};

#endif
