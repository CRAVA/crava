#include "rplib/distributionsco2.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/demmodelling.h"

DistributionsCO2::DistributionsCO2(const DistributionWithTrend       * distr_temperature,
                                   const DistributionWithTrend       * distr_pore_pressure,
                                   std::vector<double>               & alpha)
: DistributionsFluid()
{
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
  alpha_               = alpha;
}

DistributionsCO2::DistributionsCO2(const DistributionsCO2 & dist)
: DistributionsFluid(dist)
{
  distr_temperature_   = dist.distr_temperature_->Clone();
  distr_pore_pressure_ = dist.distr_pore_pressure_->Clone();
  alpha_               = dist.alpha_;
}

DistributionsCO2::~DistributionsCO2()
{
}

DistributionsFluid *
DistributionsCO2::Clone() const
{
  return new DistributionsCO2(*this);
}

Fluid *
DistributionsCO2::GenerateSample(const std::vector<double> & trend_params) const
{
  double  temperature   = distr_temperature_->ReSample(trend_params[0], trend_params[1]);
  double  pore_pressure = distr_pore_pressure_->ReSample(trend_params[0], trend_params[1]);
  Fluid * fluid         = new CO2(temperature, pore_pressure);

  return  fluid;
}

bool
DistributionsCO2::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsCO2::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Fluid *
DistributionsCO2::UpdateSample(double                      corr_param,
                               bool                        param_is_time,
                               const std::vector<double> & /*trend*/,
                               const Fluid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Fluid * updated_sample = NULL; // = GetSample();

  return updated_sample;
}
