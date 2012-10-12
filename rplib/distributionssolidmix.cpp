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
  bool has_distribution = false;

  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_solids; i++) {

    if(distr_solid_[i]->HasDistribution() == true)
      has_distribution = true;

    else if(distr_vol_frac_[i] != NULL && distr_vol_frac_[i]->GetIsDistribution() == true)
      has_distribution = true;
  }

  return has_distribution;
}

std::vector<bool>
DistributionsSolidMix::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_solids; i++) {
    std::vector<bool> solid_trend  = distr_solid_[i]->HasTrend();

    std::vector<bool> volume_trend(2,false);

    if(distr_vol_frac_[i] != NULL)
       volume_trend = distr_vol_frac_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {

      if(solid_trend[j] == true)
        has_trend[j] = true;

      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  return has_trend;
}

Solid *
DistributionsSolidMix::UpdateSample(const std::vector< double > &/*corr*/,
                                    const Solid                 & /*solid*/) const {

  return NULL;
}
