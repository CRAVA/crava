#include "rplib/distributionsdryrockmix.h"

#include "rplib/dryrockmix.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/demmodelling.h"

#include <cassert>

#include "src/definitions.h"

#include "nrlib/random/distribution.hpp"

DistributionsDryRockMix::DistributionsDryRockMix(std::vector< DistributionsDryRock * >        & distr_dryrock,
                                                 std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                                 DEMTools::MixMethod                            mix_method,
                                                 std::vector<double>                          & alpha)
: DistributionsDryRock()
{
  assert(distr_dryrock.size() == distr_vol_frac.size());
  distr_dryrock_      = distr_dryrock;
  distr_vol_frac_   = distr_vol_frac;
  mix_method_       = mix_method;
  alpha_            = alpha;
}

DistributionsDryRockMix::DistributionsDryRockMix(const DistributionsDryRockMix & dist)
: DistributionsDryRock(dist),
  mix_method_(dist.mix_method_)
{
  for(size_t i=0; i<dist.distr_dryrock_.size(); i++)
    distr_dryrock_.push_back(dist.distr_dryrock_[i]);

  for(size_t i=0; i<dist.distr_vol_frac_.size(); i++) {
    if(dist.distr_vol_frac_[i] != NULL)
      distr_vol_frac_.push_back(dist.distr_vol_frac_[i]->Clone());
    else
      distr_vol_frac_.push_back(NULL);
  }

  alpha_ = dist.alpha_;
}

DistributionsDryRockMix::~DistributionsDryRockMix()
{
  for(size_t i=0; i<distr_dryrock_.size(); i++)
    delete distr_dryrock_[i];

  for(size_t i=0; i<distr_vol_frac_.size(); i++) {
    if(distr_vol_frac_[i] != NULL && distr_vol_frac_[i]->GetIsShared() == false)
      delete distr_vol_frac_[i];
  }
}

DistributionsDryRock *
DistributionsDryRockMix::Clone() const
{
  return new DistributionsDryRockMix(*this);
}

DryRock *
DistributionsDryRockMix::GenerateSample(const std::vector<double> & trend_params) const
{
  size_t n_dryrocks = distr_dryrock_.size();

  std::vector<double> u(n_dryrocks, RMISSING);
  for(size_t i=0; i<n_dryrocks; i++) {
    if(distr_vol_frac_[i] != NULL)
      u[i] = NRLib::Random::Unif01();
  }

  std::vector<DryRock *> dryrock_samples(n_dryrocks);

  for(size_t i = 0; i < n_dryrocks; ++i)
    dryrock_samples[i] = distr_dryrock_[i]->GenerateSample(trend_params);

  DryRock * dryrock_mixed = GetSample(u, trend_params, dryrock_samples);

  // Deep copy taken by constructor of DryRockMix, hence delete dryrock here:
  for(size_t i = 0; i < n_dryrocks; ++i)
    delete dryrock_samples[i];

  return dryrock_mixed;
}

DryRock *
DistributionsDryRockMix::GetSample(const std::vector<double>    & u,
                                   const std::vector<double>    & trend_params,
                                   const std::vector<DryRock *> & dryrock_samples) const
{
  size_t n_dryrocks = dryrock_samples.size();

  std::vector<double> volume_fraction(n_dryrocks, 0.0);

  size_t missing_index = n_dryrocks;

  for(size_t i = 0; i < n_dryrocks; ++i) {

    if (u[i] != RMISSING)
      volume_fraction[i] = distr_vol_frac_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    else
      missing_index    = i;
  }

  if (missing_index != n_dryrocks) {
    double sum = 0.0;

    for (size_t i = 0; i < n_dryrocks; ++i)
      sum += volume_fraction[i];

    volume_fraction[missing_index] = 1.0 - sum;
  }

  DryRock * dryrock_mixed = new DryRockMix(dryrock_samples, volume_fraction, u, mix_method_);

  return dryrock_mixed;
}

bool
DistributionsDryRockMix::HasDistribution() const
{
  bool has_distribution = false;

  size_t n_dryrocks = distr_dryrock_.size();

  for(size_t i=0; i<n_dryrocks; i++) {

    if(distr_dryrock_[i]->HasDistribution() == true) {
      has_distribution = true;
      break;
    }

    else if(distr_vol_frac_[i] != NULL && distr_vol_frac_[i]->GetIsDistribution() == true) {
      has_distribution = true;
      break;
    }
  }

  return has_distribution;
}

std::vector<bool>
DistributionsDryRockMix::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  size_t n_dryrocks = distr_dryrock_.size();

  for(size_t i=0; i<n_dryrocks; i++) {
    std::vector<bool> dryrock_trend  = distr_dryrock_[i]->HasTrend();

    std::vector<bool> volume_trend(2,false);

    if(distr_vol_frac_[i] != NULL)
       volume_trend = distr_vol_frac_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {

      if(dryrock_trend[j] == true)
        has_trend[j] = true;

      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  return has_trend;
}

DryRock *
DistributionsDryRockMix::UpdateSample(double                      corr_param,
                                      bool                        param_is_time,
                                      const std::vector<double> & trend,
                                      const DryRock             * sample) const {

  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);

  const DryRockMix * core_sample = dynamic_cast<const DryRockMix *>(sample);

  std::vector<DryRock *> updated_sub_dryrocks(distr_dryrock_.size());
  for(size_t i = 0; i<distr_dryrock_.size(); i++){
     updated_sub_dryrocks[i] = distr_dryrock_[i]->UpdateSample(corr_param,
                                                               param_is_time,
                                                               trend,
                                                               core_sample->GetSubDryRock(i));
  }

  DryRock * updated_sample = GetSample(u, trend, updated_sub_dryrocks);

  for(size_t i = 0; i<distr_dryrock_.size(); i++)
    delete updated_sub_dryrocks[i];

  return updated_sample;
}
