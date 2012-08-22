#ifndef RPLIB_DISTRIBUTIONSBRINE_H
#define RPLIB_DISTRIBUTIONSBRINE_H

#include "rplib/distributionsfluid.h"
#include "rplib/brine.h"
#include "rplib/distributionsbrineevolution.h"

#include "nrlib/random/distribution.hpp"

class DistributionsBrine : public DistributionsFluid {
public:

  DistributionsBrine(NRLib::Distribution<double> * distr_temperature,
                     NRLib::Distribution<double> * distr_pore_pressure,
                     NRLib::Distribution<double> * distr_salinity,
                     DistributionsBrineEvolution * distr_evolution = NULL);

  virtual ~DistributionsBrine();

  virtual Fluid * GenerateSample(const std::vector<double> & /*trend_params*/) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

private:
  NRLib::Distribution<double> * distr_salinity_;      // Pointer to external object.
  NRLib::Distribution<double> * distr_temperature_;   // Pointer to external object.
  NRLib::Distribution<double> * distr_pore_pressure_; // Pointer to external object.
  DistributionsBrineEvolution * distr_evolution_;     // Pointer to external object.

};

#endif
