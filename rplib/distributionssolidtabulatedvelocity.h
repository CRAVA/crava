#ifndef RPLIB_DISTRIBUTIONS_SOLID_TABULATED_VELOCITY_H
#define RPLIB_DISTRIBUTIONS_SOLID_TABULATED_VELOCITY_H

#include "rplib/solid.h"
#include "rplib/distributionssolid.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/tabulated.h"


// Abstract class for holding all t = 0 distribution functions for solid parameters.
// One derived class for each solid model, the latter specified in a parallel, derived Solid class.
// The class must be able to produce an object of the specific Solid class.
class DistributionsSolidTabulatedVelocity : public DistributionsSolid {
public:

  //NB: Class is not completed
  DistributionsSolidTabulatedVelocity(const DistributionWithTrend * vp,
                                      const DistributionWithTrend * vs,
                                      const DistributionWithTrend * density,
                                      double                        corr_vp_vs,
                                      double                        corr_vp_density,
                                      double                        corr_vs_density);

  virtual                       ~DistributionsSolidTabulatedVelocity();

  // Solid is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Solid               * GenerateSample(const std::vector<double> & trend_params) const;

  virtual bool                  HasDistribution() const;

  virtual std::vector<bool>     HasTrend() const;

protected:
  virtual Solid *               UpdateSample(const std::vector< double > & /*corr*/,
                                             const Solid                 & /*solid*/) const;

private:
  Solid                       * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend   * vp_;
  const DistributionWithTrend   * vs_;
  const DistributionWithTrend   * density_;
  double                          corr_vp_vs_;
  double                          corr_vp_density_;
  double                          corr_vs_density_;
  Tabulated                     * tabulated_;
  bool                            has_distribution_;
  std::vector<bool>               has_trend_;
};

#endif
