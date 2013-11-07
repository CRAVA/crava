#include "rplib/distributionsdryrockmix.h"

#include "rplib/distributionssolid.h"

#include "rplib/dryrockmix.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/demmodelling.h"

#include <cassert>

#include "src/definitions.h"

#include "nrlib/random/distribution.hpp"

//This file contains two classes DistributionsDryRockMix and DistributionsDryRockMixOfDryRockAndSolid.

//-------------------------------------- DistributionsDryRockMix ---------------------------------------------------------

DistributionsDryRockMix::DistributionsDryRockMix(std::vector< DistributionsDryRock * >        & distr_dryrock,
                                                 std::vector< DistributionWithTrend * >       & distr_vol_frac,
                                                 DEMTools::MixMethod                            mix_method,
                                                 std::vector<double>                          & alpha)
: DistributionsDryRock(),
  distr_dryrock_(distr_dryrock.size(), NULL),
  distr_vol_frac_(distr_vol_frac.size(), NULL),
  mix_method_(mix_method)
{
  assert(distr_dryrock.size() == distr_vol_frac.size());

  for (size_t i = 0; i < distr_dryrock.size(); ++i)
    distr_dryrock_[i] = distr_dryrock[i]->Clone();

  for (size_t i = 0; i < distr_vol_frac.size(); ++i) {
    if(distr_vol_frac_[i] != NULL) {
      if(distr_vol_frac[i]->GetIsShared() == false)
        distr_vol_frac_[i] = distr_vol_frac[i]->Clone();
      else
        distr_vol_frac_[i] = distr_vol_frac[i];
    }
  }
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

  for (size_t i = 0; i < distr_vol_frac_.size(); i++) {
    if (distr_vol_frac_[i] != NULL) {
      if (distr_vol_frac_[i]->GetIsShared() == false)
      delete distr_vol_frac_[i];
    }
  }
}

DistributionsDryRock *
DistributionsDryRockMix::Clone() const
{
  return new DistributionsDryRockMix(*this);
}

DryRock *
DistributionsDryRockMix::GenerateSample(const std::vector<double> & trend_params)
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
                                   const std::vector<DryRock *> & dryrock_samples)
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
                                      const DryRock             * sample)
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

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

//-------------------------------------- DistributionsDryRockMixOfDryRockAndSolid ---------------------------------------------------------


DistributionsDryRockMixOfDryRockAndSolid::DistributionsDryRockMixOfDryRockAndSolid(std::vector< DistributionsDryRock * >        & distr_dryrock,
                                                                                   std::vector< DistributionsSolid* >           & distr_solid,
                                                                                   std::vector< DistributionWithTrend * >       & distr_vol_frac_dryrock,
                                                                                   std::vector< DistributionWithTrend * >       & distr_vol_frac_solid,
                                                                                   DEMTools::MixMethod                            mix_method,
                                                                                   std::vector<double>                          & alpha)
: DistributionsDryRock(),
  mix_method_(mix_method)
{
  assert((distr_solid.size() == distr_vol_frac_solid.size()) && (distr_dryrock.size() ==  distr_vol_frac_dryrock.size()));

  size_t solid_size = distr_solid.size();
  distr_solid_.resize(solid_size);
  for(size_t i=0; i<solid_size; i++)
    distr_solid_[i] = distr_solid[i]->Clone();

  size_t dryrock_size = distr_dryrock.size();
  distr_dryrock_.resize(dryrock_size);
  for(size_t i=0; i<dryrock_size; i++)
    distr_dryrock_[i] = distr_dryrock[i]->Clone();

  distr_vol_frac_solid_.resize(solid_size,NULL);
  for(size_t i=0; i<solid_size; i++) {
    if(distr_vol_frac_solid[i] != NULL) {
      if(distr_vol_frac_solid[i]->GetIsShared() == false)
        distr_vol_frac_solid_[i] = distr_vol_frac_solid[i]->Clone();
      else
        distr_vol_frac_solid_[i] = distr_vol_frac_solid[i];
    }
  }

  distr_vol_frac_dryrock_.resize(dryrock_size,NULL);
  for(size_t i=0; i<dryrock_size; i++) {
    if(distr_vol_frac_dryrock[i] != NULL) {
      if(distr_vol_frac_dryrock[i]->GetIsShared() == false)
        distr_vol_frac_dryrock_[i] = distr_vol_frac_dryrock[i]->Clone();
      else
        distr_vol_frac_dryrock_[i] = distr_vol_frac_dryrock[i];
    }
  }

  alpha_ = alpha;               // alpha_ contains the one-year correlations for (distr_vol_fraction_solid, distr_vol_fraction_dryrock)
}

DistributionsDryRockMixOfDryRockAndSolid::DistributionsDryRockMixOfDryRockAndSolid(const DistributionsDryRockMixOfDryRockAndSolid & dist)
: DistributionsDryRock(dist),
  mix_method_(dist.mix_method_)
{
  size_t solid_size = dist.distr_solid_.size();
  size_t dryrock_size = dist.distr_dryrock_.size();

  distr_solid_.resize(solid_size);
  for(size_t i=0; i<solid_size; i++)
    distr_solid_[i] = dist.distr_solid_[i]->Clone();

  distr_dryrock_.resize(dryrock_size);
  for(size_t i=0; i<dryrock_size; i++)
    distr_dryrock_[i] = dist.distr_dryrock_[i]->Clone();

  distr_vol_frac_solid_.resize(solid_size, NULL);
  for(size_t i=0; i<solid_size; i++) {
    if(dist.distr_vol_frac_solid_[i] != NULL) {
      if(dist.distr_vol_frac_solid_[i]->GetIsShared() == false)
        distr_vol_frac_solid_[i] = dist.distr_vol_frac_solid_[i]->Clone();
      else
        distr_vol_frac_solid_[i] = dist.distr_vol_frac_solid_[i];
    }
  }

  distr_vol_frac_dryrock_.resize(dryrock_size, NULL);
  for(size_t i=0; i<dryrock_size; i++) {
    if(dist.distr_vol_frac_dryrock_[i] != NULL) {
      if(dist.distr_vol_frac_dryrock_[i]->GetIsShared() == false)
        distr_vol_frac_dryrock_[i] = dist.distr_vol_frac_dryrock_[i]->Clone();
      else
        distr_vol_frac_dryrock_[i] = dist.distr_vol_frac_dryrock_[i];
    }
  }

  alpha_       = dist.alpha_;
}

DistributionsDryRockMixOfDryRockAndSolid::~DistributionsDryRockMixOfDryRockAndSolid()
{
  for(size_t i=0; i<distr_solid_.size(); i++)
    delete distr_solid_[i];

  for(size_t i=0; i<distr_dryrock_.size(); i++)
    delete distr_dryrock_[i];

  for(size_t i=0; i<distr_vol_frac_solid_.size(); i++) {
    if(distr_vol_frac_solid_[i] != NULL) {
      if(distr_vol_frac_solid_[i]->GetIsShared() == false)
        delete distr_vol_frac_solid_[i];
    }
  }

  for(size_t i=0; i<distr_vol_frac_dryrock_.size(); i++) {
    if(distr_vol_frac_dryrock_[i] != NULL) {
      if(distr_vol_frac_dryrock_[i]->GetIsShared() == false)
        delete distr_vol_frac_dryrock_[i];
    }
  }
}

DistributionsDryRock *
DistributionsDryRockMixOfDryRockAndSolid::Clone() const
{
  return new DistributionsDryRockMixOfDryRockAndSolid(*this);
}

DryRock *
DistributionsDryRockMixOfDryRockAndSolid::GetSample(const std::vector<double>     & u,
                                                    const std::vector<double>     & trend_params,
                                                    const std::vector<Solid *>    & solid_sample,
                                                    const std::vector<DryRock *>  & dryrock_sample)
{

  size_t n_dryrocks    =   dryrock_sample.size();
  size_t n_solids      =   solid_sample.size();

  std::vector<double> volume_fraction_solid(n_solids, 0.0);

  size_t missing_index = n_solids + n_dryrocks;

  for(size_t i = 0; i < n_solids; i++) {
    if (u[i] != RMISSING)
      volume_fraction_solid[i] = distr_vol_frac_solid_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    else
      missing_index = i;
  }

  std::vector<double> volume_fraction_dryrock(n_dryrocks, 0.0);

  for(size_t i = 0; i < n_dryrocks; i++) {
    if (u[i + n_solids] != RMISSING)
      volume_fraction_dryrock[i] = distr_vol_frac_dryrock_[i]->GetQuantileValue(u[i + n_solids], trend_params[0], trend_params[1]);
    else
      missing_index = i+n_solids;
  }

  if (missing_index != n_dryrocks + n_solids) {
    double sum = 0.0;

    for (size_t i = 0; i < volume_fraction_solid.size(); ++i)
      sum += volume_fraction_solid[i];

    for (size_t i = 0; i < volume_fraction_dryrock.size(); ++i)
      sum += volume_fraction_dryrock[i];

   if(missing_index < n_solids)
     volume_fraction_solid[missing_index] = 1.0 - sum;
   else
    volume_fraction_dryrock[missing_index - n_solids] = 1.0 - sum;
  }

  DryRock * dryrock_mixed = new DryRockMixOfDryRockAndSolid(dryrock_sample, solid_sample, volume_fraction_dryrock, volume_fraction_solid, u, mix_method_);

  return dryrock_mixed;
}

DryRock *
DistributionsDryRockMixOfDryRockAndSolid::GenerateSample(const std::vector<double> & trend_params)
{
  size_t n_dryrocks    =   distr_dryrock_.size();
  size_t n_solids      =   distr_solid_.size();

  std::vector<double> u(n_solids + n_dryrocks, RMISSING);

  for(size_t i=0; i<n_solids; i++) {
    if(distr_vol_frac_solid_[i] != NULL)
      u[i] = NRLib::Random::Unif01();
  }

  for(size_t i=0; i<n_dryrocks; i++) {
    if(distr_vol_frac_dryrock_[i] != NULL)
      u[i + n_solids] = NRLib::Random::Unif01();
  }

  std::vector<Solid *> solid_sample(n_solids);

  for(size_t i = 0; i < n_solids; ++i)
    solid_sample[i] = distr_solid_[i]->GenerateSample(trend_params);

  std::vector<DryRock *> dryrock_sample(n_dryrocks);

  for(size_t i = 0; i < n_dryrocks; ++i)
    dryrock_sample[i] = distr_dryrock_[i]->GenerateSample(trend_params);

  DryRock * dryrock_mixed = GetSample(u, trend_params, solid_sample, dryrock_sample);

  // Deep copy taken by constructor of DistributionsDryRockMixOfDryRockAndSolid, hence delete solid and dryrock here:
  for(size_t i = 0; i < n_solids; ++i)
    delete solid_sample[i];

  for(size_t i = 0; i < n_dryrocks; ++i)
    delete dryrock_sample[i];

  return dryrock_mixed;
}

bool
DistributionsDryRockMixOfDryRockAndSolid::HasDistribution() const
{
  bool has_distribution = false;

  size_t n_dryrocks = distr_dryrock_.size();
  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_dryrocks; i++) {

    if(distr_dryrock_[i]->HasDistribution() == true)
      has_distribution = true;

    else if(distr_vol_frac_dryrock_[i] != NULL && distr_vol_frac_dryrock_[i]->GetIsDistribution() == true)
      has_distribution = true;
  }
  for(size_t i=0; i<n_solids; i++) {

    if(distr_solid_[i]->HasDistribution() == true)
      has_distribution = true;

    else if(distr_vol_frac_solid_[i] != NULL && distr_vol_frac_solid_[i]->GetIsDistribution() == true)
      has_distribution = true;
  }

  return has_distribution;

}

std::vector<bool>
DistributionsDryRockMixOfDryRockAndSolid::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  size_t n_dryrocks = distr_dryrock_.size();

  for(size_t i=0; i<n_dryrocks; i++) {
    std::vector<bool> dryrock_trend  = distr_dryrock_[i]->HasTrend();

    std::vector<bool> volume_trend(2, false);

    if(distr_vol_frac_dryrock_[i] != NULL)
       volume_trend = distr_vol_frac_dryrock_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {

      if(dryrock_trend[j] == true)
        has_trend[j] = true;

      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  size_t n_solids = distr_solid_.size();

  for(size_t i=0; i<n_solids; i++) {
    std::vector<bool> solid_trend  = distr_solid_[i]->HasTrend();

    std::vector<bool> volume_trend(2,false);

    if(distr_vol_frac_solid_[i] != NULL)
       volume_trend = distr_vol_frac_solid_[i]->GetUseTrendCube();

    for(int j=0; j<2; j++) {

      if(solid_trend[j] == true)
        has_trend[j] = true;

      else if(volume_trend[j] == true)
        has_trend[j] = true;
    }
  }

  return has_trend;
}

DryRock *
DistributionsDryRockMixOfDryRockAndSolid::UpdateSample(double                      corr_param,
                                                       bool                        param_is_time,
                                                       const std::vector<double> & trend,
                                                       const DryRock             * sample)
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

  assert(typeid(*sample) == typeid(DryRockMixOfDryRockAndSolid));
  const DryRockMixOfDryRockAndSolid * core_sample = dynamic_cast<const DryRockMixOfDryRockAndSolid *>(sample);

  std::vector<Solid *>   updated_solids(distr_solid_.size());
  std::vector<DryRock *> updated_dryrocks(distr_dryrock_.size());

  for(size_t i = 0; i < distr_solid_.size(); i++) {
    updated_solids[i] = distr_solid_[i]->UpdateSample(corr_param,
                                                      param_is_time,
                                                      trend,
                                                      core_sample->GetSubSolid(i));
  }

  for(size_t i = 0; i < distr_dryrock_.size(); i++) {
    updated_dryrocks[i] = distr_dryrock_[i]->UpdateSample(corr_param,
                                                          param_is_time,
                                                          trend,
                                                          core_sample->GetSubDryRock(i));
  }

  DryRock * updated_sample = GetSample(u, trend, updated_solids, updated_dryrocks);

  for(size_t i = 0; i < distr_solid_.size(); i++)
    delete updated_solids[i];

  for(size_t i = 0; i < distr_dryrock_.size(); i++)
    delete updated_dryrocks[i];

  return updated_sample;
}

