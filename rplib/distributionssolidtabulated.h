#ifndef RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H
#define RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H

#include "rplib/solid.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"


// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolidTabulated : public DistributionsSolid {
public:

  //NB: Class is not completed
  DistributionsSolidTabulated(const DistributionWithTrend * elastic1,
                              const DistributionWithTrend * elastic2,
                              const DistributionWithTrend * density,
                              double                        corr_elastic1_elastic2,
                              double                        corr_elastic1_density,
                              double                        corr_elastic2_density,
                              DEMTools::TabulatedMethod     method,
                              std::vector<double>         & alpha);

  virtual ~DistributionsSolidTabulated();

  // Solid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Solid               * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

protected:
  virtual Solid *               UpdateSample(double                      corr_param,
                                             bool                        param_is_time,
                                             const std::vector<double> & trend,
                                             const Solid               * sample) const;

private:
  Solid                       * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend   * elastic1_;
  const DistributionWithTrend   * elastic2_;
  const DistributionWithTrend   * density_;
  double                          corr_elastic1_elastic2_;
  double                          corr_elastic1_density_;
  double                          corr_elastic2_density_;
  Tabulated                     * tabulated_;
  DEMTools::TabulatedMethod       tabulated_method_;
};

#endif
