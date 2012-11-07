#ifndef RPLIB_DISTRIBUTIONS_FLUID_TABULATED_H
#define RPLIB_DISTRIBUTIONS_FLUID_TABULATED_H

#include "rplib/fluid.h"
#include "rplib/distributionsfluid.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"

// Abstract class for holding all t = 0 distribution functions for fluid parameters.
// One derived class for each fluid model, the latter specified in a parallel, derived Fluid class.
// The class must be able to produce an object of the specific Fluid class.
class DistributionsFluidTabulated : public DistributionsFluid {
public:

  DistributionsFluidTabulated(const DistributionWithTrend * elastic,
                              const DistributionWithTrend * density,
                              double                        corr_elastic_density,
                              DEMTools::TabulatedMethod     method,
                              std::vector<double>         & alpha);

  virtual ~DistributionsFluidTabulated();

  // Fluid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Fluid               * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

  virtual Fluid               * UpdateSample(double                      corr_param,
                                             bool                        param_is_time,
                                             const std::vector<double> & trend,
                                             const Fluid               * sample) const;

private:

  Fluid                       * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend * elastic_;
  const DistributionWithTrend * density_;
  double                        corr_elastic_density_;
  Tabulated                   * tabulated_;
  DEMTools::TabulatedMethod     tabulated_method_;

};

#endif
