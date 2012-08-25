#ifndef RPLIB_DISTRIBUTIONSCO2_H
#define RPLIB_DISTRIBUTIONSCO2_H

#include "rplib/distributionsfluid.h"
#include "rplib/co2.h"

class DistributionWithTrend;
class DistributionsCO2Evolution;

class DistributionsCO2 : public DistributionsFluid {
public:

  DistributionsCO2(DistributionWithTrend       * distr_temperature,
                   DistributionWithTrend       * distr_pore_pressure,
                   DistributionsCO2Evolution   * distr_evolution = NULL);

  virtual ~DistributionsCO2();

  virtual Fluid * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

private:
  DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.
  DistributionsCO2Evolution   * distr_evolution_;     // Pointer to external object.
};

#endif
