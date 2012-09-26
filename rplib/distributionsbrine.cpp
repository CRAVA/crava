#include "rplib/distributionsbrine.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/distributionsbrineevolution.h"

DistributionsBrine::DistributionsBrine(const DistributionWithTrend        * distr_temperature,
                                       const DistributionWithTrend        * distr_pore_pressure,
                                       const DistributionWithTrend        * distr_salinity,
                                             DistributionsBrineEvolution  * distr_evolution)
: DistributionsFluid()
{
  distr_salinity_      = distr_salinity;
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
  distr_evolution_     = distr_evolution;
}

DistributionsBrine::~DistributionsBrine(){}

Fluid *
DistributionsBrine::GenerateSample(const std::vector<double> & trend_params) const
{
  double  salinity      = distr_salinity_->ReSample(trend_params[0], trend_params[1]);
  double  temperature   = distr_temperature_->ReSample(trend_params[0], trend_params[1]);
  double  pore_pressure = distr_pore_pressure_->ReSample(trend_params[0], trend_params[1]);
  Fluid * fluid         = new Brine(salinity, temperature, pore_pressure, distr_evolution_);

  return  fluid;
}

bool
DistributionsBrine::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsBrine::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Fluid *
DistributionsBrine::UpdateSample(const std::vector< double > & /*corr*/,
                                 const Fluid                 & /*fluid*/) const {

  return NULL;
}
