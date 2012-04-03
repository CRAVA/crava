#ifndef SOLID_H
#define SOLID_H

#include <vector>
#include "rplib/distributionssolidevolution.h"

// Abstract solid class.
// Each derived class has parallel classes derived from DistributionsSolid and DistributionsSolidEvolve.
class Solid {
public:

  Solid() {}
  virtual ~Solid() {}

  virtual void ComputeElasticParams(double & k, double & mu, double & rho) const = 0;

  // Solid is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Derived class Evolve implementation should always start with casting and assert,
  // analogous to what is done for Fluid.
  virtual Solid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Solid * >       & solid,
                         const DistributionsSolidEvolution  * dist_solid_evolve) const = 0;

protected:

};

#endif
