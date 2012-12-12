#include "rplib/distributionsrockdem.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/rockdem.h"
#include "rplib/demmodelling.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/statistics/statistics.hpp"
#include <nrlib/flens/nrlib_flens.hpp>

#include "src/definitions.h"

#include <cassert>
#include <vector>


DistributionsRockDEM::DistributionsRockDEM(DistributionsSolid                           * distr_solid,
                                           std::vector< DistributionsFluid *>           & distr_fluid,
                                           std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                           std::vector< DistributionWithTrend * >       & distr_incl_concentration,
                                           const std::vector<double>                    & alpha,
                                           const std::vector<double>                    & s_min,
                                           const std::vector<double>                    & s_max)
: DistributionsRock()
{
  assert( distr_incl_spectrum.size() + 1 == distr_incl_concentration.size() );
  assert( distr_incl_spectrum.size()     == distr_fluid.size());

  int n_inclusions = static_cast<int>(distr_incl_spectrum.size());

  distr_solid_ = distr_solid->Clone();

  distr_fluid_.resize(n_inclusions);
  for(int i=0; i<n_inclusions; i++)
    distr_fluid_[i] = distr_fluid[i]->Clone();

  distr_incl_spectrum_.resize(n_inclusions);
  for(int i=0; i<n_inclusions; i++) {
    if(distr_incl_spectrum[i]->GetIsShared() == false)
      distr_incl_spectrum_[i] = distr_incl_spectrum[i]->Clone();
    else
      distr_incl_spectrum_[i] = distr_incl_spectrum[i];
  }

  distr_incl_concentration_.resize(n_inclusions+1, NULL);
  for(int i=0; i<n_inclusions+1; i++) {
    if(distr_incl_concentration[i] != NULL) {
      if(distr_incl_concentration[i]->GetIsShared() == false)
        distr_incl_concentration_[i] = distr_incl_concentration[i]->Clone();
      else
        distr_incl_concentration_[i] = distr_incl_concentration[i];
    }
  }

  alpha_                    = alpha;               // alpha_ contains the one-year correlations for (inclusion_spectrums, inclusion_concentrations)
  s_min_                    = s_min;
  s_max_                    = s_max;

  SampleVpVsRhoExpectationAndCovariance(expectation_old_, covariance_old_);

  SetupExpectationAndCovariances(expectation_,
                                 covariance_,
                                 tabulated_s0_,
                                 tabulated_s1_,
                                 s_min_,
                                 s_max_);
}

DistributionsRockDEM::DistributionsRockDEM(const DistributionsRockDEM & dist)
: DistributionsRock(dist),
  expectation_old_(dist.expectation_old_),
  covariance_old_(dist.covariance_old_)
{
  int n_inclusions = static_cast<int>(dist.distr_incl_spectrum_.size());

  distr_solid_ = dist.distr_solid_->Clone();

  distr_fluid_.resize(n_inclusions);
  for(int i=0; i<n_inclusions; i++)
    distr_fluid_[i] = dist.distr_fluid_[i]->Clone();

  distr_incl_spectrum_.resize(n_inclusions);
  for(int i=0; i<n_inclusions; i++) {
    if(dist.distr_incl_spectrum_[i]->GetIsShared() == false)
      distr_incl_spectrum_[i] = dist.distr_incl_spectrum_[i]->Clone();
    else
      distr_incl_spectrum_[i] = dist.distr_incl_spectrum_[i];
  }

  distr_incl_concentration_.resize(n_inclusions+1, NULL);
  for(int i=0; i<n_inclusions+1; i++) {
    if(dist.distr_incl_concentration_[i] != NULL) {
      if(dist.distr_incl_concentration_[i]->GetIsShared() == false)
        distr_incl_concentration_[i] = dist.distr_incl_concentration_[i]->Clone();
      else
        distr_incl_concentration_[i] = dist.distr_incl_concentration_[i];
    }
  }

  alpha_       = dist.alpha_;
  s_min_       = dist.s_min_;
  s_max_       = dist.s_max_;
  expectation_ = dist.expectation_;
  covariance_  = dist.covariance_;
}

DistributionsRockDEM::~DistributionsRockDEM()
{
  delete distr_solid_;

  for(size_t i=0; i<distr_fluid_.size(); i++)
    delete distr_fluid_[i];

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

DistributionsRock *
DistributionsRockDEM::Clone() const
{
  return new DistributionsRockDEM(*this);
}


Rock *
DistributionsRockDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  Solid * solid = distr_solid_->GenerateSample(trend_params);

  std::vector< Fluid* > fluid(distr_fluid_.size());
  for (size_t i = 0; i < fluid.size(); ++i)
    fluid[i] = distr_fluid_[i]->GenerateSample(trend_params);

  size_t n_incl = distr_incl_spectrum_.size();

  std::vector<double> u(n_incl+n_incl+1, RMISSING);
  for(size_t i=0; i<n_incl; i++) {
    if (distr_incl_concentration_[i] != NULL)
      u[i + n_incl] = NRLib::Random::Unif01();

    u[i] = NRLib::Random::Unif01();
  }

  //last element incl check
  if (distr_incl_concentration_.back() != NULL)
    u.back() = NRLib::Random::Unif01();

  Rock * new_rock = GetSample(u, trend_params, solid, fluid);

  // Deep copy taken by constructor of RockInclusion, hence delete
  // solid and fluid here:
  delete solid;
  for (size_t i = 0; i < fluid.size(); ++i)
    delete fluid[i];

  return new_rock;
}

bool
DistributionsRockDEM::HasDistribution() const
{
  if (distr_solid_->HasDistribution())
      return true;

  for (size_t i = 0; i < distr_fluid_.size(); ++i) {
    if (distr_fluid_[i]->HasDistribution())
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
DistributionsRockDEM::HasTrend() const
{

  std::vector<bool> has_trend(2, false);

  std::vector<bool> solid_trend     = distr_solid_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    const std::vector<bool>& incl_trend         = distr_incl_spectrum_[i]->GetUseTrendCube();
    std::vector<bool> incl_conc(2, false);
    if (distr_incl_concentration_[i] != NULL)
      incl_conc  = distr_incl_concentration_[i]->GetUseTrendCube();
    const std::vector<bool>& fluid_trend_inc    = distr_fluid_[i]->HasTrend();

    for(size_t j = 0; j < 2; ++j) {
      if (solid_trend[j] || fluid_trend_inc[j] || incl_trend[j] || incl_conc[j])
        has_trend[j] = true;
    }
  }

  //check last element in inclusion
  std::vector<bool> incl_conc(2, false);
  if (distr_incl_concentration_.back() != NULL)
    incl_conc  = distr_incl_concentration_.back()->GetUseTrendCube();

  for(size_t j = 0; j < 2; ++j) {
    if (incl_conc[j])
      has_trend[j] = true;
  }

  return has_trend;
}

void
DistributionsRockDEM::SampleVpVsRhoExpectationAndCovariance(std::vector<double>   & expectation,
                                                            NRLib::Grid2D<double> & covariance)
{
  int nsamples = 100;

  std::vector< NRLib::Vector > m(3);
  NRLib::Vector vp(nsamples);
  NRLib::Vector vs(nsamples);
  NRLib::Vector rho(nsamples);

  std::vector<double> dummy(2, 0.0);
  for (int i = 0; i < nsamples; ++i) {
    Rock * rock = GenerateSample(dummy);
    rock->GetSeismicParams(vp(i), vs(i), rho(i));
    delete rock;
  }
  m[0] = vp;
  m[1] = vs;
  m[2] = rho;

  expectation.resize(3, 0.0);
  covariance.Resize(3,3,0.0);
  for (int i = 0; i < 3; ++i) {
    expectation[i] = NRLib::Mean(m[i]);
    for (int j = 0; j < 3; ++j)
      covariance(i,j) = NRLib::Cov(m[i], m[j]);
  }
}

Rock *
DistributionsRockDEM::GetSample(const std::vector<double>  & u,
                                const std::vector<double>  & trend_params,
                                const Solid                * solid,
                                const std::vector< Fluid *>& fluid) const
{
  size_t n_incl = distr_incl_spectrum_.size();
  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl + 1, 0.0);

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

  return new RockDEM(solid, fluid, inclusion_spectrum, inclusion_concentration, u);
}

Rock *
DistributionsRockDEM::UpdateSample(double                      corr_param,
                                   bool                        param_is_time,
                                   const std::vector<double> & trend,
                                   const Rock                * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, alpha_);

  assert(typeid(*sample) == typeid(RockDEM));
  const RockDEM * core_sample = dynamic_cast<const RockDEM *>(sample);

  Solid * updated_solid = distr_solid_->UpdateSample(corr_param,
                                                     param_is_time,
                                                     trend,
                                                     core_sample->GetSolid());

  std::vector<Fluid *> updated_fluid(distr_fluid_.size());
  for (size_t i = 0; i < updated_fluid.size(); ++i) {
  updated_fluid[i]      = distr_fluid_[i]->UpdateSample(corr_param,
                                                        param_is_time,
                                                        trend,
                                                        core_sample->GetFluid(i));
  }

  Rock * updated_sample = GetSample(u, trend, updated_solid, updated_fluid);

  return updated_sample;
}
