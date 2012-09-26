#ifndef RPLIB_DISTRIBUTIONS_FLUID_TABULATED_VELOCITY_H
#define RPLIB_DISTRIBUTIONS_FLUID_TABULATED_VELOCITY_H

#include "rplib/fluid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrend.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluidTabulatedVelocity : public DistributionsFluid {
public:

  //NB: Class is not completed
  DistributionsFluidTabulatedVelocity(const DistributionWithTrend * vp,
                                      const DistributionWithTrend * density,
                                      double                        corr_vp_density);

  virtual                       ~DistributionsFluidTabulatedVelocity();

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid               * GenerateSample(const std::vector<double> & /*trend_params*/) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

protected:
  virtual Fluid               * UpdateSample(const std::vector< double > & /*corr*/,
                                             const Fluid                 & /*fluid*/) const;

private:
  const DistributionWithTrend * vp_;
  const DistributionWithTrend * density_;
  double                        corr_vp_density_;
  bool                          has_distribution_;
};

#endif
