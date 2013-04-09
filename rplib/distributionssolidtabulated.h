#ifndef RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H
#define RPLIB_DISTRIBUTIONS_SOLID_TABULATED_H

#include "rplib/distributionssolid.h"
#include "rplib/demmodelling.h"

class Solid;
class DistributionWithTrend;
class Tabulated;


// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolidTabulated : public DistributionsSolid {
public:

  //NB: Class is not completed
  DistributionsSolidTabulated(DistributionWithTrend       * elastic1,
                              DistributionWithTrend       * elastic2,
                              DistributionWithTrend       * density,
                              double                        corr_elastic1_elastic2,
                              double                        corr_elastic1_density,
                              double                        corr_elastic2_density,
                              DEMTools::TabulatedMethod     method,
                              std::vector<double>         & alpha);

  DistributionsSolidTabulated(const DistributionsSolidTabulated & dist);

  virtual ~DistributionsSolidTabulated();

  // Solid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.

  virtual DistributionsSolid  * Clone() const;

  virtual Solid               * GenerateSample(const std::vector<double> & trend_params);

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

protected:
  virtual Solid *               UpdateSample(double                      corr_param,
                                             bool                        param_is_time,
                                             const std::vector<double> & trend,
                                             const Solid               * sample);

private:
  Solid                       * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params);

  DistributionWithTrend         * elastic1_;
  DistributionWithTrend         * elastic2_;
  DistributionWithTrend         * density_;
  double                          corr_elastic1_elastic2_;
  double                          corr_elastic1_density_;
  double                          corr_elastic2_density_;
  Tabulated                     * tabulated_;
  DEMTools::TabulatedMethod       tabulated_method_;
};

#endif
