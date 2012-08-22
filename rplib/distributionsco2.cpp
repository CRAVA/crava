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
DistributionsCO2::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  double temperature   = distr_temperature_->Draw();
  double pore_pressure = distr_pore_pressure_->Draw();
  Fluid * fluid = new CO2(temperature, pore_pressure, distr_evolution_);

  return fluid;
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
