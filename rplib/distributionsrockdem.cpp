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
                                           DistributionsFluid                           * distr_fluid,
                                           std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                           std::vector< DistributionWithTrend * >       & distr_incl_concentration)
: DistributionsRock()
{
  assert( distr_incl_spectrum.size() == distr_incl_concentration.size() );

  distr_solid_              = distr_solid;
  distr_fluid_              = distr_fluid;
  distr_incl_spectrum_      = distr_incl_spectrum;
  distr_incl_concentration_ = distr_incl_concentration;

  SampleVpVsRhoExpectationAndCovariance(expectation_old_, covariance_old_);
}

DistributionsRockDEM::~DistributionsRockDEM(){}

Rock *
DistributionsRockDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  Solid * solid = distr_solid_->GenerateSample(trend_params);
  Fluid * fluid = distr_fluid_->GenerateSample(trend_params);
  size_t n_incl = distr_incl_spectrum_.size();

  std::vector<double> u(n_incl+n_incl);
  for(size_t i=0; i<n_incl+n_incl; i++)
    u[i] = NRLib::Random::Unif01();

  Rock * new_rock = GetSample(u, trend_params, solid, fluid);

  // Deep copy taken by constructor of RockInclusion, hence delete
  // solid and fluid here:
  delete solid;
  delete fluid;

  return new_rock;
}

Pdf3D *
DistributionsRockDEM::GeneratePdf(void) const
{
  Pdf3D * pdf3D = NULL;
  return pdf3D;
}

bool
DistributionsRockDEM::HasDistribution() const
{
  if (distr_solid_->HasDistribution() || distr_fluid_->HasDistribution())
      return true;

  // loop over inclusion and spectrum
  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    if (distr_incl_spectrum_[i]->GetIsDistribution() || distr_incl_concentration_[i]->GetIsDistribution())
      return true;
  }

  return false;
}

std::vector<bool>
DistributionsRockDEM::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> solid_trend     = distr_solid_->HasTrend();
  std::vector<bool> fluid_trend     = distr_fluid_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    std::vector<bool> incl_trend = distr_incl_spectrum_[i]->GetUseTrendCube();
    std::vector<bool> incl_conc  = distr_incl_concentration_[i]->GetUseTrendCube();

    for(size_t j = 0; j < 2; ++j) {
      if (solid_trend[j] || fluid_trend[j] || incl_trend[j] || incl_conc[j])
        has_trend[j] = true;
    }
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
                                const Fluid                * fluid) const
{
  size_t n_incl = distr_incl_spectrum_.size();
  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl);

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i] = distr_incl_spectrum_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    inclusion_concentration[i] = distr_incl_concentration_[i]->GetQuantileValue(u[i + n_incl], trend_params[0], trend_params[1]);
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
  DEMTools::UpdateU(u, corr_param, param_is_time);

  assert(typeid(sample) == typeid(RockDEM));
  const RockDEM * core_sample = dynamic_cast<const RockDEM *>(sample);

  Solid * updated_solid = distr_solid_->UpdateSample(corr_param,
                                                     param_is_time,
                                                     trend,
                                                     core_sample->GetSolid());
  Fluid * updated_fluid = distr_fluid_->UpdateSample(corr_param,
                                                     param_is_time,
                                                     trend,
                                                     core_sample->GetFluid());

  Rock * updated_sample = GetSample(u, trend, updated_solid, updated_fluid);

  return updated_sample;
}
