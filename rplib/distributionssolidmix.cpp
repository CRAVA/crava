#include "rplib/distributionssolidmix.h"

#include <cassert>

DistributionsSolidMix::DistributionsSolidMix(std::vector< DistributionsSolid * >            distr_solid,
                                             std::vector< NRLib::Distribution<double> * >   distr_vol_frac,
                                             DEMTools::MixMethod                            mix_method,
                                             DistributionsSolidMixEvolution               * distr_evolution)
: DistributionsSolid()
{
  assert(distr_solid.size() == distr_vol_frac.size());
  distr_solid_      = distr_solid;
  distr_vol_frac_   = distr_vol_frac;
  mix_method_       = mix_method;
  distr_evolution_  = distr_evolution;
}

DistributionsSolidMix::~DistributionsSolidMix(){}

Solid *
DistributionsSolidMix::GenerateSample() const
{

  size_t n_solids = distr_solid_.size();
  std::vector<Solid*> solid(n_solids);
  std::vector<double> volume_fraction(n_solids);

  for(size_t i = 0; i < n_solids; ++i) {
    solid[i] = distr_solid_[i]->GenerateSample();
    volume_fraction[i] = distr_vol_frac_[i]->Draw();
  }
  Solid * solid_mixed = new SolidMixed(solid, volume_fraction, mix_method_, distr_evolution_);

  // Deep copy taken by constructor of SolidMixed, hence delete solid here:
  for(size_t i = 0; i < n_solids; ++i)
    delete solid[i];

  return solid_mixed;
}
