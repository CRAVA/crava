#ifndef RPLIB_DISTRIBUTIONSBRINE_H
#define RPLIB_DISTRIBUTIONSBRINE_H

#include "rplib/distributionsfluid.h"
#include "rplib/brine.h"
#include "rplib/distributionsbrineevolution.h"

class DistributionsBrine : public DistributionsFluid {
public:

  DistributionsBrine(NRLib::Distribution<double> * distr_temperature,
                     NRLib::Distribution<double> * distr_pore_pressure,
                     NRLib::Distribution<double> * distr_salinity,
                     DistributionsBrineEvolution * distr_evolution = NULL)
  : DistributionsFluid()
  {
    distr_salinity_      = distr_salinity;
    distr_temperature_   = distr_temperature;
    distr_pore_pressure_ = distr_pore_pressure;
    distr_evolution_     = distr_evolution;
  }

  virtual ~DistributionsBrine(){}

  virtual Fluid * GenerateSample() const {
    double salinity      = distr_salinity_->Draw();
    double temperature   = distr_temperature_->Draw();
    double pore_pressure = distr_pore_pressure_->Draw();
    Fluid * fluid = new Brine(salinity, temperature, pore_pressure, distr_evolution_);

    return fluid;
  }


private:
  NRLib::Distribution<double> * distr_salinity_;
  NRLib::Distribution<double> * distr_temperature_;
  NRLib::Distribution<double> * distr_pore_pressure_;
  DistributionsBrineEvolution * distr_evolution_;
};

#endif
