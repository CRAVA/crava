#ifndef RPLIB_DISTRIBUTIONSSOLIDMIX_H
#define RPLIB_DISTRIBUTIONSSOLIDMIX_H

#include "rplib/distributionssolid.h"
#include "rplib/solidmixed.h"
#include "rplib/distributionssolidmixevolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsSolidMix : public DistributionsSolid {
public:

  DistributionsSolidMix(std::vector< DistributionsSolid * >            distr_solid,
                        std::vector< NRLib::Distribution<double> * >   distr_vol_frac,
                        DistributionsSolidMixEvolution               * distr_evolution = NULL);

  virtual ~DistributionsSolidMix();

  virtual Solid * GenerateSample() const;

private:
  std::vector< DistributionsSolid * >            distr_solid_;     // Pointers to external objects.
  std::vector< NRLib::Distribution<double> * >   distr_vol_frac_;  // Pointers to external objects.
  DistributionsSolidMixEvolution               * distr_evolution_; // Pointer to external object.
};

#endif
