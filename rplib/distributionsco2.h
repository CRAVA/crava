#ifndef RPLIB_DISTRIBUTIONSCO2_H
#define RPLIB_DISTRIBUTIONSCO2_H

#include "rplib/distributionsfluid.h"
#include "rplib/co2.h"

class DistributionWithTrend;

class DistributionsCO2 : public DistributionsFluid {
public:

                                      DistributionsCO2(const DistributionWithTrend       * distr_temperature,
                                                     const DistributionWithTrend       * distr_pore_pressure);

  virtual                             ~DistributionsCO2();

  virtual Fluid                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                        HasDistribution() const;

  virtual std::vector<bool>           HasTrend() const;

  virtual Fluid *                     UpdateSample(double                      corr_param,
                                                   bool                        param_is_time,
                                                   const std::vector<double> & trend,
                                                   const Fluid               * sample) const;
protected:

private:
  const DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  const DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.
};

#endif
