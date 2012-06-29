#include "rplib/distributionsbrine.h"

DistributionsBrine::DistributionsBrine(NRLib::Distribution<double> * distr_temperature,
                                       NRLib::Distribution<double> * distr_pore_pressure,
                                       NRLib::Distribution<double> * distr_salinity,
                                       DistributionsBrineEvolution * distr_evolution)
: DistributionsFluid()
{
  distr_salinity_      = distr_salinity;
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
  distr_evolution_     = distr_evolution;
}

DistributionsBrine::~DistributionsBrine(){}

Fluid *
DistributionsBrine::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  double salinity      = distr_salinity_->Draw();
  double temperature   = distr_temperature_->Draw();
  double pore_pressure = distr_pore_pressure_->Draw();
  Fluid * fluid = new Brine(salinity, temperature, pore_pressure, distr_evolution_);

  return fluid;
}
