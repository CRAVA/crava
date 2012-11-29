#ifndef RPLIB_DISTRIBUTIONS_ROCK_TABULATED_H
#define RPLIB_DISTRIBUTIONS_ROCK_TABULATED_H

#include "rplib/distributionsrock.h"
#include "rplib/demmodelling.h"

class Rock;
class Tabulated;
class DistributionWithTrend;

class DistributionsRockTabulated : public DistributionsRock {
public:

  DistributionsRockTabulated(const DistributionWithTrend * elastic1,
                             const DistributionWithTrend * elastic2,
                             const DistributionWithTrend * density,
                             double                        corr_elastic1_elastic2,
                             double                        corr_elastic1_density,
                             double                        corr_elastic2_density,
                             DEMTools::TabulatedMethod     method,
                             const std::vector<double>   & alpha,
                             const std::vector<double>   & s_min,
                             const std::vector<double>   & s_max);

  DistributionsRockTabulated(const DistributionsRockTabulated & dist);

  virtual ~DistributionsRockTabulated();

  virtual DistributionsRock        * Clone() const;

  // Rock is an abstract class, hence pointer must be used here. Allocated memory (using new) MUST be deleted by caller.
  virtual Rock                     * GenerateSample(const std::vector<double> & trend_params) const;

  virtual Rock                     * UpdateSample(double                      corr_param,
                                                  bool                        param_is_time,
                                                  const std::vector<double> & trend,
                                                  const Rock                * sample)       const;

  virtual bool                       HasDistribution() const;

  virtual std::vector<bool>          HasTrend() const;

  virtual bool                       GetIsOkForBounding() const;

private:

  Rock                             * GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const;

  const DistributionWithTrend * elastic1_;
  const DistributionWithTrend * elastic2_;
  const DistributionWithTrend * density_;
  double                        corr_elastic1_elastic2_;
  double                        corr_elastic1_density_;
  double                        corr_elastic2_density_;
  Tabulated                   * tabulated_;
  DEMTools::TabulatedMethod     tabulated_method_;
};

#endif
