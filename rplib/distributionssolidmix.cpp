#include "rplib/distributionssolidmix.h"

#include "rplib/solidmix.h"
#include "rplib/distributionwithtrend.h"

#include <cassert>

DistributionsSolidMix::DistributionsSolidMix(std::vector< DistributionsSolid * >          & distr_solid,
                                             std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                             DEMTools::MixMethod                            mix_method)
: DistributionsSolid()
{
  assert(distr_solid.size() == distr_vol_frac.size());
  distr_solid_      = distr_solid;
  distr_vol_frac_   = distr_vol_frac;
  mix_method_       = mix_method;
}

DistributionsSolidMix::~DistributionsSolidMix(){}

Solid *
DistributionsSolidMix::GenerateSample(const std::vector<double> & trend_params) const
{

  size_t n_solids = distr_solid_.size();
  std::vector<Solid*> solid(n_solids);
  std::vector<double> volume_fraction(n_solids, 0.0);

  size_t missing_index = n_solids;
  for(size_t i = 0; i < n_solids; ++i) {
    solid[i] = distr_solid_[i]->GenerateSample(trend_params);
    if (distr_vol_frac_[i])
      volume_fraction[i] = distr_vol_frac_[i]->ReSample(trend_params[0], trend_params[1]);
    else
      missing_index    = i;
  }

  if (missing_index != n_solids) {
    double sum = 0.0;
    for (size_t i = 0; i < volume_fraction.size(); ++i)
      sum += volume_fraction[i];

    volume_fraction[missing_index] = 1.0 - sum;
  }

  Solid * solid_mixed = new SolidMix(solid, volume_fraction, mix_method_);

  // Deep copy taken by constructor of SolidMix, hence delete solid here:
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

Solid *
DistributionsSolidMix::UpdateSample(const std::vector< double > &/*corr*/,
                                    const Solid                 & /*solid*/) const {

  return NULL;
}
