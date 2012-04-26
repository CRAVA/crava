#include "rplib/distributionsco2.h"

DistributionsCO2::DistributionsCO2(NRLib::Distribution<double> * distr_temperature,
                                   NRLib::Distribution<double> * distr_pore_pressure,
                                   DistributionsCO2Evolution   * distr_evolution)
: DistributionsFluid()
{
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
  distr_evolution_     = distr_evolution;
}

DistributionsCO2::~DistributionsCO2(){}

Fluid *
DistributionsCO2::GenerateSample() const
{
  double temperature   = distr_temperature_->Draw();
  double pore_pressure = distr_pore_pressure_->Draw();
  Fluid * fluid = new CO2(temperature, pore_pressure, distr_evolution_);

  return fluid;
}
