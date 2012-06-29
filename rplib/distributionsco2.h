#ifndef RPLIB_DISTRIBUTIONSCO2_H
#define RPLIB_DISTRIBUTIONSCO2_H

#include "rplib/distributionsfluid.h"
#include "rplib/co2.h"
#include "rplib/distributionsco2evolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsCO2 : public DistributionsFluid {
public:

  DistributionsCO2(NRLib::Distribution<double> * distr_temperature,
                   NRLib::Distribution<double> * distr_pore_pressure,
                   DistributionsCO2Evolution   * distr_evolution = NULL);

  virtual ~DistributionsCO2();

  virtual Fluid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

private:
  NRLib::Distribution<double> * distr_temperature_;   // Pointer to external object.
  NRLib::Distribution<double> * distr_pore_pressure_; // Pointer to external object.
  DistributionsCO2Evolution   * distr_evolution_;     // Pointer to external object.
};

#endif
