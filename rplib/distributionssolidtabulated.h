#ifndef RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H
#define RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H

#include "rplib/solid.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionwithtrend.h"

// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolidTabulated : public DistributionsSolid {
public:

  //NB: Class is not completed
  DistributionsSolidTabulated(const DistributionWithTrend * vp,
                              const DistributionWithTrend * vs,
                              const DistributionWithTrend * density,
                              const DistributionWithTrend * corr_vp_vs,
                              const DistributionWithTrend * corr_vp_density,
                              const DistributionWithTrend * corr_vs_density);

  virtual ~DistributionsSolidTabulated();

  // Solid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Solid * GenerateSample() const;

private:
  const DistributionWithTrend * vp_;
  const DistributionWithTrend * vs_;
  const DistributionWithTrend * density_;
  const DistributionWithTrend * corr_vp_vs_;
  const DistributionWithTrend * corr_vp_density_;
  const DistributionWithTrend * corr_vs_density_;
};

#endif
