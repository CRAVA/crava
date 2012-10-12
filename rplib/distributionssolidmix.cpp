#include "rplib/distributionssolidmix.h"

#include "rplib/solidmix.h"
#include "rplib/distributionwithtrend.h"

#include <cassert>

#include "src/definitions.h"

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

  std::vector<double> u(n_solids, RMISSING);
  for(size_t i=0; i<n_solids; i++) {
    if(distr_vol_frac_[i] != NULL)
      u[i] = NRLib::Random::Unif01();
  }

  std::vector<Solid *> solid_samples(n_solids);

  for(size_t i = 0; i < n_solids; ++i)
    solid_samples[i] = distr_solid_[i]->GenerateSample(trend_params);

  Solid * solid_mixed = GetSample(u, trend_params, solid_samples);

  // Deep copy taken by constructor of SolidMix, hence delete solid here:
  for(size_t i = 0; i < n_solids; ++i)
    delete solid_samples[i];

  return solid_mixed;
}

Solid *
DistributionsSolidMix::GetSample(const std::vector<double>  & u,
                                 const std::vector<double>  & trend_params,
                                 const std::vector<Solid *> & solid_samples) const
{

  size_t n_solids = solid_samples.size();

  std::vector<double> volume_fraction(n_solids, 0.0);

  size_t missing_index = n_solids;

  for(size_t i = 0; i < n_solids; ++i) {

    if (u[i] == RMISSING)
      volume_fraction[i] = distr_vol_frac_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    else
      missing_index    = i;
  }

  if (missing_index != n_solids) {
    double sum = 0.0;

    for (size_t i = 0; i < n_solids; ++i)
      sum += volume_fraction[i];

    volume_fraction[missing_index] = 1.0 - sum;
  }

  Solid * solid_mixed = new SolidMix(solid_samples, volume_fraction, u, mix_method_);

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
