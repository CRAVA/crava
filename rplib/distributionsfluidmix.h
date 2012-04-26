#ifndef RPLIB_DISTRIBUTIONSFLUIDMIX_H
#define RPLIB_DISTRIBUTIONSFLUIDMIX_H

#include "rplib/distributionsfluid.h"
#include "rplib/fluidmixed.h"
#include "rplib/distributionsfluidmixevolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsFluidMix : public DistributionsFluid {
public:

  DistributionsFluidMix(std::vector< DistributionsFluid * >            distr_fluid,
                        std::vector< NRLib::Distribution<double> * >   distr_vol_frac,
                        DistributionsFluidMixEvolution               * distr_evolution = NULL);

  virtual ~DistributionsFluidMix();

  virtual Fluid * GenerateSample() const;

private:
  std::vector< DistributionsFluid * >            distr_fluid_;     // Pointers to external objects.
  std::vector< NRLib::Distribution<double> * >   distr_vol_frac_;  // Pointers to external objects.
  DistributionsFluidMixEvolution               * distr_evolution_; // Pointer to external object.
};

#endif
