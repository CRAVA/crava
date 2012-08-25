#include "rplib/distributionssolidmix.h"

#include "rplib/distributionssolidmixevolution.h"
#include "rplib/solidmixed.h"
#include "rplib/distributionwithtrend.h"

#include <cassert>

DistributionsSolidMix::DistributionsSolidMix(std::vector< DistributionsSolid * >          & distr_solid,
                                             std::vector< DistributionWithTrend * >       & distr_vol_frac,
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
DistributionsSolidMix::GenerateSample(const std::vector<double> & trend_params) const
{

  size_t n_solids = distr_solid_.size();
  std::vector<Solid*> solid(n_solids);
  std::vector<double> volume_fraction(n_solids);

  for(size_t i = 0; i < n_solids; ++i) {
    solid[i] = distr_solid_[i]->GenerateSample(trend_params);
    volume_fraction[i] = distr_vol_frac_[i]->ReSample(trend_params[0], trend_params[1]);
  }
  Solid * solid_mixed = new SolidMixed(solid, volume_fraction, mix_method_, distr_evolution_);

  // Deep copy taken by constructor of SolidMixed, hence delete solid here:
  for(size_t i = 0; i < n_solids; ++i)
    delete solid[i];

  return solid_mixed;
}

bool
DistributionsSolidMix::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsSolidMix::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}
