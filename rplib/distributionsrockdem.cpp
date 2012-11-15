#include "rplib/distributionsrockdem.h"

#include "rplib/distributionwithtrend.h"

#include "nrlib/statistics/statistics.hpp"

#include "rplib/distributionssolid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/rockdem.h"
#include "rplib/demmodelling.h"

#include "nrlib/random/distribution.hpp"

#include <cassert>

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

  distr_solid_              = distr_solid;
  distr_fluid_              = distr_fluid;
  distr_incl_spectrum_      = distr_incl_spectrum;
  distr_incl_concentration_ = distr_incl_concentration;
  alpha_                    = alpha;
  s_min_                    = s_min;
  s_max_                    = s_max;

  SampleVpVsRhoExpectationAndCovariance(expectation_old_, covariance_old_);
}

DistributionsRockDEM::DistributionsRockDEM(const DistributionsRockDEM & dist)
: DistributionsRock(dist),
  expectation_old_(dist.expectation_old_),
  covariance_old_(dist.covariance_old_)
{
  distr_solid_ = dist.distr_solid_->Clone();

  for(size_t i=0; i<dist.distr_fluid_.size(); i++)
    distr_fluid_.push_back(dist.distr_fluid_[i]);

  for(size_t i=0; i<dist.distr_incl_spectrum_.size(); i++)
    distr_incl_spectrum_.push_back(dist.distr_incl_spectrum_[i]->Clone());

  for(size_t i=0; i<dist.distr_incl_concentration_.size(); i++)
    distr_incl_concentration_.push_back(dist.distr_incl_concentration_[i]->Clone());

  alpha_ = dist.alpha_;
  s_min_ = dist.s_min_;
  s_max_ = dist.s_max_;
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
    if(distr_incl_concentration_[i]->GetIsShared() == false)
      delete distr_incl_concentration_[i];
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

  std::vector<double> u(n_incl+n_incl+1);
  for(size_t i=0; i<n_incl+n_incl+1; i++)
    u[i] = NRLib::Random::Unif01();

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

  // loop over inclusion and spectrum
  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    if (distr_incl_spectrum_[i]->GetIsDistribution() || distr_incl_concentration_[i]->GetIsDistribution())
      return true;
  }

  // check last element
  if (distr_incl_concentration_.back()->GetIsDistribution())
    return true;

  return false;
}

std::vector<bool>
DistributionsRockDEM::HasTrend() const
{

  std::vector<bool> has_trend(2, false);

  std::vector<bool> solid_trend     = distr_solid_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    const std::vector<bool>& incl_trend         = distr_incl_spectrum_[i]->GetUseTrendCube();
    const std::vector<bool>& incl_conc          = distr_incl_concentration_[i]->GetUseTrendCube();
    const std::vector<bool>& fluid_trend_inc    = distr_fluid_[i]->HasTrend();

    for(size_t j = 0; j < 2; ++j) {
      if (solid_trend[j] || fluid_trend_inc[j] || incl_trend[j] || incl_conc[j])
        has_trend[j] = true;
    }
  }

  //check last element
  const std::vector<bool>& incl_conc  = distr_incl_concentration_.back()->GetUseTrendCube();
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
  std::vector<double> inclusion_concentration(n_incl + 1);

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i] = distr_incl_spectrum_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    inclusion_concentration[i] = distr_incl_concentration_[i]->GetQuantileValue(u[i + n_incl], trend_params[0], trend_params[1]);
  }

  inclusion_concentration.back() = distr_incl_concentration_.back()->GetQuantileValue(u.back(), trend_params[0], trend_params[1]);

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

  assert(typeid(sample) == typeid(RockDEM));
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
