#ifndef RPLIB_DISTRIBUTIONS_FLUID_TABULATED_H
#define RPLIB_DISTRIBUTIONS_FLUID_TABULATED_H

#include "rplib/fluid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrend.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluidTabulated : public DistributionsFluid {
public:

  //NB: Class is not completed
  DistributionsFluidTabulated(const DistributionWithTrend * vp,
                              const DistributionWithTrend * density,
                              const DistributionWithTrend * corr_vp_density);

  virtual ~DistributionsFluidTabulated();

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid * GenerateSample() const;

private:
  const DistributionWithTrend * vp_;
  const DistributionWithTrend * density_;
  const DistributionWithTrend * corr_vp_density_;
};

#endif
