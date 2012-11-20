#ifndef RPLIB_DISTRIBUTIONS_FLUID_CO2_H
#define RPLIB_DISTRIBUTIONS_FLUID_CO2_H

#include "rplib/distributionsfluid.h"

class DistributionWithTrend;
class FluidCO2;

class DistributionsFluidCO2 : public DistributionsFluid {
public:

  DistributionsFluidCO2(const DistributionWithTrend       * distr_temperature,
                        const DistributionWithTrend       * distr_pore_pressure,
                        std::vector<double>               & alpha);

  DistributionsFluidCO2(const DistributionsFluidCO2 & dist);

  virtual ~DistributionsFluidCO2();

  virtual DistributionsFluid        * Clone() const;

  virtual Fluid                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                        HasDistribution() const;

  virtual std::vector<bool>           HasTrend() const;

  virtual Fluid *                     UpdateSample(double                      corr_param,
                                                   bool                        param_is_time,
                                                   const std::vector<double> & trend,
                                                   const Fluid               * sample) const;

private:
  Fluid                             * GetSample(const std::vector<double> & u,
                                                const std::vector<double> & trend_params) const;

  const DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  const DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.
};

#endif
