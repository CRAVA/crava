#ifndef RPLIB_FLUIDMIXED_H
#define RPLIB_FLUIDMIXED_H

#include "rplib/fluid.h"
#include "rplib/distributionsfluidmixevolution.h"
#include "rplib/demmodelling.h"

#include "nrlib/exception/exception.hpp"

#include <cassert>
#include <vector>
#include <numeric>

class FluidMixed : public Fluid {
public:

  FluidMixed(const std::vector<Fluid*>      & fluid,
             const std::vector<double>      & volume_fraction,
             DistributionsFluidMixEvolution * distr_evolution = NULL);

  virtual ~FluidMixed();

  // Assignment operator, not yet implemented.
  /*FluidMixed& operator=(const FluidMixed& rhs);*/

  virtual Fluid * Clone() const;

  virtual void ComputeElasticParams(double /*temp*/, double /*pore_pressure*/);

  virtual void GetElasticParams(double& k, double& rho) const;

  virtual Fluid * Evolve(const std::vector<int>       & delta_time,
                         const std::vector< Fluid * > & fluid) const;

private:
  //Copy constructor for getting base class variables , used by Clone:
  FluidMixed(const FluidMixed & rhs) : Fluid(rhs) {}

  std::vector<Fluid*> fluid_;                        // Owned and deleted by this class.
  std::vector<double> volume_fraction_;
  DistributionsFluidMixEvolution * distr_evolution_; // Pointer to external object.

  double k_;
  double rho_;

};

#endif
