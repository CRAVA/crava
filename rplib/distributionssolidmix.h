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
                        DistributionsSolidMixEvolution               * distr_evolution = NULL)
  : DistributionsSolid()
  {
    assert(distr_solid.size() == distr_vol_frac.size());
    distr_solid_      = distr_solid;
    distr_vol_frac_   = distr_vol_frac;
    distr_evolution_  = distr_evolution;
  }

  virtual ~DistributionsSolidMix(){}

  virtual Solid * GenerateSample() const {

    size_t n_solids = distr_solid_.size();
    std::vector<Solid*> solid(n_solids);
    std::vector<double> volume_fraction(n_solids);

    for(size_t i = 0; i < n_solids; ++i) {
      solid[i] = distr_solid_[i]->GenerateSample();
      volume_fraction[i] = distr_vol_frac_[i]->Draw();
    }
    Solid * solid_mixed = new SolidMixed(solid, volume_fraction, distr_evolution_);

    // Deep copy taken by constructor of SolidMixed, hence delete solid here:
    for(size_t i = 0; i < n_solids; ++i)
      delete solid[i];

    return solid_mixed;
  }

private:
  std::vector< DistributionsSolid * >            distr_solid_;
  std::vector< NRLib::Distribution<double> * >   distr_vol_frac_;
  DistributionsSolidMixEvolution               * distr_evolution_;
};

#endif
