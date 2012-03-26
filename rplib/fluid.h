#ifndef FLUID_H
#define FLUID_H

#include <string>
#include <vector>
#include "rplib/distributionsfluidevolution.h"

// Abstract fluid class.
// Each derived class has parallel classes derived from DistributionsFluid and DistributionsFluidEvolve.
class Fluid {
public:

  Fluid();
  virtual ~Fluid();

  virtual void ComputeElasticParams(const double   temp,
                                    const double   pore_pressure,
                                    double       & k,
                                    double       & rho) const = 0;

  // Fluid is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Derived class Evolve implementation should always start with casting and assert.
  // Example from derived class Brine:
          //const DistributionsBrineEvolution * dist_brine_evolve =
          //      dynamic_cast<const DistributionsBrineEvolution*>(dist_fluid_evolve);
          //assert(dist_brine_evolve != NULL);
          //assert(delta_time.size() == fluid.size() + 1);
  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Fluid * >       & fluid,
                         const DistributionsFluidEvolution  * dist_fluid_evolve) const = 0;


protected:
};

#endif
