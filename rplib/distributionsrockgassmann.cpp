#include "rplib/distributionsrockgassmann.h"
#include "rplib/distributionsdryrock.h"
#include "rplib/distributionsfluid.h"
#include "rplib/rockgassmann.h"

#include <typeinfo>
#include <cassert>
#include <typeinfo>

DistributionsRockGassmann::DistributionsRockGassmann(const DistributionsDryRock            * distr_dryrock,
                                                     const DistributionsFluid              * distr_fluid,
                                                     const std::vector<double>             & s_min,
                                                     const std::vector<double>             & s_max)
: DistributionsRock()
{

  distr_dryrock_ = distr_dryrock->Clone();
  distr_fluid_   = distr_fluid->Clone();

  s_min_                    = s_min;
  s_max_                    = s_max;
}

DistributionsRockGassmann::DistributionsRockGassmann(const DistributionsRockGassmann & dist)
: DistributionsRock(dist)
{

  distr_dryrock_ = dist.distr_dryrock_->Clone();
  distr_fluid_   = dist.distr_fluid_->Clone();

  alpha_       = dist.alpha_;
  s_min_       = dist.s_min_;
  s_max_       = dist.s_max_;
  expectation_ = dist.expectation_;
  covariance_  = dist.covariance_;
}

DistributionsRockGassmann::~DistributionsRockGassmann()
{
  delete distr_fluid_;
  delete distr_dryrock_;
}

DistributionsRock *
DistributionsRockGassmann::Clone() const
{
  return new DistributionsRockGassmann(*this);
}


Rock *
DistributionsRockGassmann::GenerateSamplePrivate(const std::vector<double> & trend_params)
{
  DryRock * dryrock = distr_dryrock_->GenerateSample(trend_params);
  Fluid   * fluid   = distr_fluid_->GenerateSample(trend_params);

  Rock * new_rock = GetSample(dryrock, fluid);

  // Deep copy taken by constructor of RockInclusion, hence delete
  // dryrock and fluid here:
  delete fluid;
  delete dryrock;

  return new_rock;
}

bool
DistributionsRockGassmann::HasDistribution() const
{
  if (distr_dryrock_->HasDistribution() || distr_fluid_->HasDistribution())
      return true;

  return false;
}

std::vector<bool>
DistributionsRockGassmann::HasTrend() const
{

  std::vector<bool> has_trend(2, false);

  std::vector<bool> dryrock_trend     = distr_dryrock_->HasTrend();
  std::vector<bool> fluid_trend       = distr_fluid_->HasTrend();

  for(size_t j = 0; j < 2; ++j) {
    if (dryrock_trend[j] || fluid_trend[j])
      has_trend[j] = true;
  }
  return has_trend;
}

Rock *
DistributionsRockGassmann::GetSample(const DryRock              * dryrock,
                                     const Fluid                * fluid)
{

  return new RockGassmann(fluid, dryrock);
}

Rock *
DistributionsRockGassmann::UpdateSample(double                      corr_param,
                                        bool                        param_is_time,
                                        const std::vector<double> & trend,
                                        const Rock                * sample)
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

  assert(typeid(*sample) == typeid(RockGassmann));
  const RockGassmann * core_sample = dynamic_cast<const RockGassmann *>(sample);

  DryRock * updated_dryrock = distr_dryrock_->UpdateSample(corr_param,
                                                           param_is_time,
                                                           trend,
                                                           core_sample->GetDryRock());

  Fluid * updated_fluid = distr_fluid_->UpdateSample(corr_param,
                                                     param_is_time,
                                                     trend,
                                                     core_sample->GetFluid());

  Rock * updated_sample = GetSample(updated_dryrock, updated_fluid);

  delete updated_fluid;
  delete updated_dryrock;

  return updated_sample;
}
