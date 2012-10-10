#ifndef RPLIB_DISTRIBUTIONSBRINE_H
#define RPLIB_DISTRIBUTIONSBRINE_H

#include "rplib/distributionsfluid.h"
#include "rplib/brine.h"

class DistributionWithTrend;

class DistributionsBrine : public DistributionsFluid {
public:

  DistributionsBrine(const DistributionWithTrend * distr_temperature,
                     const DistributionWithTrend * distr_pore_pressure,
                     const DistributionWithTrend * distr_salinity);

  virtual ~DistributionsBrine();

  virtual Fluid                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                        HasDistribution() const;

  virtual std::vector<bool>           HasTrend() const;

protected:
  virtual Fluid *                     UpdateSample(const std::vector< double > & /*corr*/,
                                                   const Fluid                 & /*fluid*/) const;

private:
  const DistributionWithTrend       * distr_salinity_;      // Pointer to external object.
  const DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  const DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.

};

#endif
