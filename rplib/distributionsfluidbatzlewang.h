#ifndef RPLIB_DISTRIBUTIONS_FLUID_BATZLE_WANG_H
#define RPLIB_DISTRIBUTIONS_FLUID_BATZLE_WANG_H

#include "rplib/distributionsfluid.h"

class Fluid;

class DistributionWithTrend;

class DistributionsFluidBatzleWang : public DistributionsFluid {
public:

  DistributionsFluidBatzleWang(const DistributionWithTrend * distr_temperature,
                               const DistributionWithTrend * distr_pore_pressure,
                               const DistributionWithTrend * distr_salinity,
                               std::vector<double>         & alpha);

  virtual ~DistributionsFluidBatzleWang();

  virtual Fluid                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                        HasDistribution() const;

  virtual std::vector<bool>           HasTrend() const;

protected:
  virtual Fluid *                     UpdateSample(double                      corr_param,
                                                   bool                        param_is_time,
                                                   const std::vector<double> & trend,
                                                   const Fluid               * sample) const;

private:

  Fluid                             * GetSample(const std::vector<double> & u,
                                                const std::vector<double> & trend_params) const;

  const DistributionWithTrend       * distr_salinity_;      // Pointer to external object.
  const DistributionWithTrend       * distr_temperature_;   // Pointer to external object.
  const DistributionWithTrend       * distr_pore_pressure_; // Pointer to external object.
};

#endif
