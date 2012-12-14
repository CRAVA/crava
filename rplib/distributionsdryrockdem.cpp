#include "rplib/distributionsdryrockdem.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/dryrockdem.h"
#include "rplib/demmodelling.h"

#include "src/definitions.h"

#include "nrlib/random/distribution.hpp"

#include <cassert>

DistributionsDryRockDEM::DistributionsDryRockDEM(DistributionsDryRock                         * distr_dryrock,
                                                 std::vector<DistributionsDryRock*>           & distr_dryrock_inc,
                                                 std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                                 std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                                                 std::vector<double>                          & alpha)
: DistributionsDryRock()
{
  assert( distr_incl_spectrum.size() + 1 == distr_incl_concentration.size() );

  distr_dryrock_              = distr_dryrock;
  distr_dryrock_inc_          = distr_dryrock_inc;
  distr_incl_spectrum_        = distr_incl_spectrum;
  distr_incl_concentration_   = distr_incl_concentration;
  alpha_                      = alpha;   // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions

}

DistributionsDryRockDEM::DistributionsDryRockDEM(const DistributionsDryRockDEM & dist)
: DistributionsDryRock(dist)
{
  distr_dryrock_ = dist.distr_dryrock_->Clone();

  for(size_t i=0; i<dist.distr_dryrock_inc_.size(); i++)
    distr_dryrock_inc_.push_back(dist.distr_dryrock_inc_[i]->Clone());

  for(size_t i=0; i<dist.distr_incl_spectrum_.size(); i++)
    distr_incl_spectrum_.push_back(dist.distr_incl_spectrum_[i]->Clone());

  for(size_t i=0; i<dist.distr_incl_concentration_.size(); i++) {
    if (dist.distr_incl_concentration_[i] != NULL)
      distr_incl_concentration_.push_back(dist.distr_incl_concentration_[i]->Clone());
    else
      distr_incl_concentration_.push_back(NULL);
  }

  alpha_ = dist.alpha_;   // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions

}

DistributionsDryRockDEM::~DistributionsDryRockDEM()
{
  delete distr_dryrock_;

  for(size_t i=0; i<distr_dryrock_inc_.size(); i++)
    delete distr_dryrock_inc_[i];

  for(size_t i=0; i<distr_incl_spectrum_.size(); i++) {
    if(distr_incl_spectrum_[i]->GetIsShared() == false)
      delete distr_incl_spectrum_[i];
  }

 for(size_t i=0; i<distr_incl_concentration_.size(); i++) {
    if (distr_incl_concentration_[i] != NULL && distr_incl_concentration_[i]->GetIsShared() == false)
      delete distr_incl_concentration_[i];
  }
}

DistributionsDryRock *
DistributionsDryRockDEM::Clone() const
{
  return new DistributionsDryRockDEM(*this);
}

DryRock *
DistributionsDryRockDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  DryRock * dryrock     = distr_dryrock_->GenerateSample(trend_params);

  std::vector< DryRock* > dryrock_inc(distr_dryrock_inc_.size());
  for (size_t i = 0; i < dryrock_inc.size(); ++i)
    dryrock_inc[i] = distr_dryrock_inc_[i]->GenerateSample(trend_params);

  size_t  n_incl    = distr_incl_spectrum_.size();

  std::vector<double> u(n_incl+n_incl+1, RMISSING);
  for(size_t i=0; i<n_incl; i++) {
    if (distr_incl_concentration_[i] != NULL)
      u[i + n_incl] = NRLib::Random::Unif01();

    u[i] = NRLib::Random::Unif01();
  }

  //last element incl check
  if (distr_incl_concentration_.back() != NULL)
    u.back() = NRLib::Random::Unif01();

  DryRock * new_dryrock = GetSample(u, trend_params, dryrock, dryrock_inc);

  // Deep copy taken by constructor of DryRockDEM, hence delete
  // dryrock and dryrock_inc here:
  delete dryrock;
  for (size_t i = 0; i < dryrock_inc.size(); ++i)
    delete dryrock_inc[i];

  return new_dryrock;
}

bool
DistributionsDryRockDEM::HasDistribution() const
{

  if (distr_dryrock_->HasDistribution())
      return true;

  for (size_t i = 0; i < distr_dryrock_inc_.size(); ++i) {
    if (distr_dryrock_inc_[i]->HasDistribution())
      return true;
  }

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    if (distr_incl_spectrum_[i]->GetIsDistribution())
      return true;
  }

  for(size_t i=0; i<distr_incl_concentration_.size(); i++) {
    if (distr_incl_concentration_[i] != NULL && distr_incl_concentration_[i]->GetIsDistribution())
      return true;
  }

  return false;
}

std::vector<bool>
DistributionsDryRockDEM::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> dryrock_trend     = distr_dryrock_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    const std::vector<bool>& incl_trend         = distr_incl_spectrum_[i]->GetUseTrendCube();

    std::vector<bool> incl_conc(2, false);
    if (distr_incl_concentration_[i] != NULL)
      incl_conc          = distr_incl_concentration_[i]->GetUseTrendCube();

    const std::vector<bool>& dryrock_trend_inc    = distr_dryrock_inc_[i]->HasTrend();

    for(size_t j = 0; j < 2; ++j) {
      if (dryrock_trend[j] || dryrock_trend_inc[j] || incl_trend[j] || incl_conc[j])
        has_trend[j] = true;
    }
  }

  std::vector<bool> incl_conc(2, false);
  if (distr_incl_concentration_.back() != NULL)
    incl_conc  = distr_incl_concentration_.back()->GetUseTrendCube();

  for(size_t j = 0; j < 2; ++j) {
    if (incl_conc[j])
      has_trend[j] = true;
  }

  return has_trend;

}

DryRock *
DistributionsDryRockDEM::UpdateSample(double                      corr_param,
                                      bool                        param_is_time,
                                      const std::vector<double> & trend,
                                      const DryRock             * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);

  assert(typeid(*sample) == typeid(DryRockDEM));
  const DryRockDEM * core_sample = dynamic_cast<const DryRockDEM *>(sample);

  DryRock * updated_dryrock_host = distr_dryrock_->UpdateSample(corr_param,
                                                          param_is_time,
                                                          trend,
                                                          core_sample->GetDryRockHost());
  std::vector<DryRock *> updated_dryrock_inc(distr_dryrock_inc_.size());
  for (size_t i = 0; i < updated_dryrock_inc.size(); ++i) {
    updated_dryrock_inc[i] = distr_dryrock_inc_[i]->UpdateSample(corr_param,
                                                             param_is_time,
                                                             trend,
                                                             core_sample->GetDryRockInclusion(i));
  }

  DryRock * updated_sample = GetSample(u, trend, updated_dryrock_host, updated_dryrock_inc);

  return updated_sample;
}

DryRock *
DistributionsDryRockDEM::GetSample(const std::vector<double>    & u,
                                   const std::vector<double>    & trend_params,
                                   const DryRock                * dryrock,
                                   const std::vector< DryRock* >& dryrock_inc) const
{
  size_t  n_incl = distr_incl_spectrum_.size();
  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl+1);

  size_t missing_index = n_incl + 1;

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i] = distr_incl_spectrum_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    if (distr_incl_concentration_[i] != NULL)
      inclusion_concentration[i] = distr_incl_concentration_[i]->GetQuantileValue(u[i + n_incl], trend_params[0], trend_params[1]);
    else
      missing_index = i;
  }

  if (distr_incl_concentration_.back() != NULL)
    inclusion_concentration.back() = distr_incl_concentration_.back()->GetQuantileValue(u.back(), trend_params[0], trend_params[1]);
  else
    missing_index = inclusion_concentration.size() - 1;

  if (missing_index != n_incl + 1) {

    double sum = 0.0;

    for (size_t i = 0; i < inclusion_concentration.size(); ++i)
      sum += inclusion_concentration[i];

    inclusion_concentration[missing_index] = 1.0 - sum;
  }

  return new DryRockDEM(dryrock, dryrock_inc, inclusion_spectrum, inclusion_concentration, u);
}
