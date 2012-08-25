#ifndef RPLIB_DISTRIBUTIONSBRINE_H
#define RPLIB_DISTRIBUTIONSBRINE_H

#include "rplib/distributionsfluid.h"
#include "rplib/brine.h"

class DistributionWithTrend;
class DistributionsBrineEvolution;

class DistributionsBrine : public DistributionsFluid {
public:

  DistributionsBrine(DistributionWithTrend * distr_temperature,
                     DistributionWithTrend * distr_pore_pressure,
                     DistributionWithTrend * distr_salinity,
                     DistributionsBrineEvolution * distr_evolution = NULL);

  virtual ~DistributionsBrine();

  virtual Fluid * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

private:
  DistributionWithTrend       * distr_salinity_;      // Pointer to external object.
  DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.
  DistributionsBrineEvolution * distr_evolution_;     // Pointer to external object.

};

#endif
