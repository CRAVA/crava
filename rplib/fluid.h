#ifndef FLUID_H
#define FLUID_H

#include <string>
#include <vector>
#include "rplib/distributionsfluidevolution.h"

// Abstract fluid class.
// Each derived class has parallel classes derived from DistributionsFluidT0 and DistributionsFluidEvolve.
class Fluid {
public:

  Fluid(double temp, double pore_pressure);
  virtual ~Fluid();

  virtual void ComputeElasticParams(double & k, double & rho) const = 0;

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Fluid * >       & fluid,
                         const DistributionsFluidEvolution  * dist_fluid_evolve) const = 0;

  void SetCommonParams(const double temp, const double pore_pressure){
    temp_          = temp;
    pore_pressure_ = pore_pressure;
  }

protected:
  // Sampled values of parameters common to all fluid classes.
  double temp_;
  double pore_pressure_;
};

#endif
