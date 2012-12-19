#include "rplib/distributionssoliddem.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/soliddem.h"
#include "rplib/demmodelling.h"

#include "src/definitions.h"

#include "nrlib/random/distribution.hpp"

#include <cassert>

DistributionsSolidDEM::DistributionsSolidDEM(DistributionsSolid                           * distr_solid,
                                             std::vector<DistributionsSolid*>             & distr_solid_inc,
                                             std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                             std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                                             std::vector<double>                          & alpha)
: DistributionsSolid(),
  distr_solid_inc_(distr_solid_inc.size(), NULL),
  distr_incl_spectrum_(distr_incl_spectrum.size(), NULL),
  distr_incl_concentration_(distr_incl_concentration.size(), NULL)

{
  assert( distr_incl_spectrum.size() + 1 == distr_incl_concentration.size() );

  distr_solid      = distr_solid->Clone();

  for (size_t i = 0; i < distr_solid_inc.size(); ++i)
    distr_solid_inc_[i]      = distr_solid_inc[i]->Clone();

  for (size_t i = 0; i < distr_incl_spectrum.size(); ++i) {
    if(distr_incl_spectrum[i]->GetIsShared() == false)
      distr_incl_spectrum_[i] = distr_incl_spectrum[i]->Clone();
    else
      distr_incl_spectrum_[i] = distr_incl_spectrum[i];
  }

  for (size_t i = 0; i < distr_incl_concentration.size(); ++i) {
    if(distr_incl_concentration[i] != NULL) {
      if(distr_incl_concentration[i]->GetIsShared() == false)
        distr_incl_concentration_[i] = distr_incl_concentration[i]->Clone();
      else
        distr_incl_concentration_[i] = distr_incl_concentration[i];
    }
  }

  alpha_                    = alpha;   // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions

}

DistributionsSolidDEM::DistributionsSolidDEM(const DistributionsSolidDEM & dist)
: DistributionsSolid(dist)
{
  distr_solid_ = dist.distr_solid_->Clone();

  size_t n_inclusions = dist.distr_solid_inc_.size();

  distr_solid_inc_.resize(n_inclusions);
  for(size_t i=0; i<n_inclusions; i++)
    distr_solid_inc_[i] = dist.distr_solid_inc_[i]->Clone();

  distr_incl_spectrum_.resize(n_inclusions);
  for(size_t i=0; i<n_inclusions; i++) {
    if(dist.distr_incl_spectrum_[i]->GetIsShared() == false)
      distr_incl_spectrum_[i] = dist.distr_incl_spectrum_[i]->Clone();
    else
      distr_incl_spectrum_[i] = dist.distr_incl_spectrum_[i];
  }

  distr_incl_concentration_.resize(n_inclusions+1, NULL);
  for(size_t i=0; i<n_inclusions+1; i++) {
    if(dist.distr_incl_concentration_[i] != NULL) {
      if(dist.distr_incl_concentration_[i]->GetIsShared() == false)
        distr_incl_concentration_[i] = dist.distr_incl_concentration_[i]->Clone();
      else
        distr_incl_concentration_[i] = dist.distr_incl_concentration_[i];
    }
  }

  alpha_ = dist.alpha_;   // Order in alpha: aspect_ratios, host_volume_fraction, inclusion_volume_fractions

}

DistributionsSolidDEM::~DistributionsSolidDEM()
{
  delete distr_solid_;

  for(size_t i=0; i<distr_solid_inc_.size(); i++)
    delete distr_solid_inc_[i];

  for(size_t i=0; i<distr_incl_spectrum_.size(); i++) {
    if(distr_incl_spectrum_[i]->GetIsShared() == false)
      delete distr_incl_spectrum_[i];
  }

  for(size_t i=0; i<distr_incl_concentration_.size(); i++) {
    if (distr_incl_concentration_[i] != NULL) {
      if(distr_incl_concentration_[i]->GetIsShared() == false)
        delete distr_incl_concentration_[i];
    }
  }
}

DistributionsSolid *
DistributionsSolidDEM::Clone() const
{
  return new DistributionsSolidDEM(*this);
}

Solid *
DistributionsSolidDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  Solid * solid     = distr_solid_->GenerateSample(trend_params);

  std::vector< Solid* > solid_inc(distr_solid_inc_.size());
  for (size_t i = 0; i < solid_inc.size(); ++i)
    solid_inc[i] = distr_solid_inc_[i]->GenerateSample(trend_params);

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

  Solid * new_solid = GetSample(u, trend_params, solid, solid_inc);

  // Deep copy taken by constructor of SolidDEM, hence delete
  // solid and solid_inc here:
  delete solid;
  for (size_t i = 0; i < solid_inc.size(); ++i)
    delete solid_inc[i];

  return new_solid;
}

bool
DistributionsSolidDEM::HasDistribution() const
{

  if (distr_solid_->HasDistribution())
      return true;

  for (size_t i = 0; i < distr_solid_inc_.size(); ++i) {
    if (distr_solid_inc_[i]->HasDistribution())
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
DistributionsSolidDEM::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> solid_trend     = distr_solid_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    const std::vector<bool>& incl_trend         = distr_incl_spectrum_[i]->GetUseTrendCube();

    std::vector<bool> incl_conc(2, false);
    if (distr_incl_concentration_[i] != NULL)
      incl_conc          = distr_incl_concentration_[i]->GetUseTrendCube();

    const std::vector<bool>& solid_trend_inc    = distr_solid_inc_[i]->HasTrend();

    for(size_t j = 0; j < 2; ++j) {
      if (solid_trend[j] || solid_trend_inc[j] || incl_trend[j] || incl_conc[j])
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

Solid *
DistributionsSolidDEM::UpdateSample(double                      corr_param,
                                    bool                        param_is_time,
                                    const std::vector<double> & trend,
                                    const Solid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);

  assert(typeid(*sample) == typeid(SolidDEM));
  const SolidDEM * core_sample = dynamic_cast<const SolidDEM *>(sample);

  Solid * updated_solid_host = distr_solid_->UpdateSample(corr_param,
                                                          param_is_time,
                                                          trend,
                                                          core_sample->GetSolidHost());
  std::vector<Solid *> updated_solid_inc(distr_solid_inc_.size());
  for (size_t i = 0; i < updated_solid_inc.size(); ++i) {
    updated_solid_inc[i] = distr_solid_inc_[i]->UpdateSample(corr_param,
                                                             param_is_time,
                                                             trend,
                                                             core_sample->GetSolidInclusion(i));
  }

  Solid * updated_sample = GetSample(u, trend, updated_solid_host, updated_solid_inc);

  delete updated_solid_host;

  for (size_t i = 0; i < updated_solid_inc.size(); ++i)
    delete updated_solid_inc[i];

  return updated_sample;
}

Solid *
DistributionsSolidDEM::GetSample(const std::vector<double>  & u,
                                 const std::vector<double>  & trend_params,
                                 const Solid                * solid,
                                 const std::vector< Solid* >& solid_inc) const
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

  return new SolidDEM(solid, solid_inc, inclusion_spectrum, inclusion_concentration, u);
}
