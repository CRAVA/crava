#include "rplib/distributionsco2.h"

#include "rplib/distributionwithtrend.h"

DistributionsCO2::DistributionsCO2(const DistributionWithTrend       * distr_temperature,
                                   const DistributionWithTrend       * distr_pore_pressure)
: DistributionsFluid()
{
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
}

DistributionsCO2::~DistributionsCO2(){}

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
DistributionsCO2::UpdateSample(const std::vector< double > & /*corr*/,
                               const Fluid                 & /*fluid*/) const {

  return NULL;
}
