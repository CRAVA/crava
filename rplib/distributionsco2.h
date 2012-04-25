#ifndef RPLIB_DISTRIBUTIONSCO2_H
#define RPLIB_DISTRIBUTIONSCO2_H

#include "rplib/distributionsfluid.h"
#include "rplib/co2.h"
#include "rplib/distributionsco2evolution.h"

class DistributionsCO2 : public DistributionsFluid {
public:

  DistributionsCO2(NRLib::Distribution<double> * distr_temperature,
                   NRLib::Distribution<double> * distr_pore_pressure,
                   DistributionsCO2Evolution   * distr_evolution = NULL)
  : DistributionsFluid()
  {
    distr_temperature_   = distr_temperature;
    distr_pore_pressure_ = distr_pore_pressure;
    distr_evolution_     = distr_evolution;
  }

  virtual ~DistributionsCO2(){}

  virtual Fluid * GenerateSample() const {
    double temperature   = distr_temperature_->Draw();
    double pore_pressure = distr_pore_pressure_->Draw();
    Fluid * fluid = new CO2(temperature, pore_pressure, distr_evolution_);

    return fluid;
  }

private:
  NRLib::Distribution<double> * distr_temperature_;
  NRLib::Distribution<double> * distr_pore_pressure_;
  DistributionsCO2Evolution   * distr_evolution_;
};

#endif
