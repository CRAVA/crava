#include "rplib/distributionssoliddem.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/soliddem.h"
#include "rplib/demmodelling.h"

#include <cassert>

DistributionsSolidDEM::DistributionsSolidDEM(DistributionsSolid                           * distr_solid,
                                             DistributionsSolid                           * distr_solid_inc,
                                             std::vector< DistributionWithTrend * >       & distr_incl_spectrum,
                                             std::vector< DistributionWithTrend * >       & distr_incl_concentration)
: DistributionsSolid()
{
  assert( distr_incl_spectrum.size() == distr_incl_concentration.size() );

  distr_solid_              = distr_solid;
  distr_solid_inc_          = distr_solid_inc;
  distr_incl_spectrum_      = distr_incl_spectrum;
  distr_incl_concentration_ = distr_incl_concentration;
}

DistributionsSolidDEM::~DistributionsSolidDEM(){}

Solid *
DistributionsSolidDEM::GenerateSample(const std::vector<double> & trend_params) const
{
  Solid * solid     = distr_solid_->GenerateSample(trend_params);
  Solid * solid_inc = distr_solid_inc_->GenerateSample(trend_params);
  size_t  n_incl    = distr_incl_spectrum_.size();

  std::vector<double> u(n_incl+n_incl);
  for(size_t i=0; i<n_incl+n_incl; i++)
    u[i] = NRLib::Random::Unif01();

  Solid * new_solid = GetSample(u, trend_params, solid, solid_inc);

  // Deep copy taken by constructor of SolidDEM, hence delete
  // solid and solid_inc here:
  delete solid;
  delete solid_inc;

  return new_solid;
}

bool
DistributionsSolidDEM::HasDistribution() const
{

  if (distr_solid_->HasDistribution() || distr_solid_inc_->HasDistribution())
      return true;

  // loop over inclusion and spectrum
  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    if (distr_incl_spectrum_[i]->GetIsDistribution() || distr_incl_concentration_[i]->GetIsDistribution())
      return true;
  }

  return false;
}

std::vector<bool>
DistributionsSolidDEM::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> solid_trend     = distr_solid_->HasTrend();
  std::vector<bool> solid_trend_inc = distr_solid_inc_->HasTrend();

  for (size_t i = 0; i < distr_incl_spectrum_.size(); ++i) {
    std::vector<bool> incl_trend = distr_incl_spectrum_[i]->GetUseTrendCube();
    std::vector<bool> incl_conc  = distr_incl_concentration_[i]->GetUseTrendCube();

    for(size_t j = 0; j < 2; ++j) {
      if (solid_trend[j] || solid_trend_inc[j] || incl_trend[j] || incl_conc[j])
        has_trend[j] = true;
    }
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

  assert(typeid(sample) == typeid(SolidDEM));
  const SolidDEM * core_sample = dynamic_cast<const SolidDEM *>(sample);

  Solid * updated_solid_host = distr_solid_->UpdateSample(corr_param,
                                                          param_is_time,
                                                          trend,
                                                          core_sample->GetSolidHost());
  Solid * updated_solid_inc = distr_solid_inc_->UpdateSample(corr_param,
                                                             param_is_time,
                                                             trend,
                                                             core_sample->GetSolidInclusion());

  Solid * updated_sample = GetSample(u, trend, updated_solid_host, updated_solid_inc);

  return updated_sample;
}

Solid *
DistributionsSolidDEM::GetSample(const std::vector<double>  & u,
                                 const std::vector<double>  & trend_params,
                                 const Solid                * solid,
                                 const Solid                * solid_inc) const
{
  size_t  n_incl = distr_incl_spectrum_.size();
  std::vector<double> inclusion_spectrum(n_incl);
  std::vector<double> inclusion_concentration(n_incl);

  for (size_t i = 0; i < n_incl; ++i) {
    inclusion_spectrum[i]      = distr_incl_spectrum_[i]->GetQuantileValue(u[i], trend_params[0], trend_params[1]);
    inclusion_concentration[i] = distr_incl_concentration_[i]->GetQuantileValue(u[i + n_incl], trend_params[0], trend_params[1]);
  }

  return new SolidDEM(solid, solid_inc, inclusion_spectrum, inclusion_concentration, u);
}
