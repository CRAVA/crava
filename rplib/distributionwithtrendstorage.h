#ifndef RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H
#define RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H

#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/trend/trendstorage.hpp"

class DistributionWithTrendStorage
{
public:
  DistributionWithTrendStorage();

  DistributionWithTrendStorage(double value,
                               bool is_sheared);

  DistributionWithTrendStorage(const NRLib::TrendStorage * trend,
                               bool                        is_sheared);

  DistributionWithTrendStorage(NRLib::Distribution<double>       * distribution,
                               const NRLib::TrendStorage         * mean,
                               const NRLib::TrendStorage         * variance,
                               bool                                is_gaussian,
                               bool                                is_sheared);

  ~DistributionWithTrendStorage();

  NRLib::TrendStorage                  * CloneMean();
  NRLib::TrendStorage                  * CloneVariance();

  const DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                       const std::vector<std::string>          & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                       std::string                             & errTxt);

  const bool                             GetIsGaussian() const;
  const bool                             GetIsDouble()   const;
  const bool                             GetIsSheared()  const;

private:
  NRLib::Distribution<double>          * distribution_;
  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
  DistributionWithTrend                * distribution_with_trend_;
  const bool                             is_gaussian_;                         // True if distribution is Gaussian
  const bool                             is_double_;                           // True if object is a double, i.e., if no distribution or trend is used
  const bool                             is_sheared_;                          // True if object is a reservoir variable that can be used for more fluids/solids/rocks/dry-rocks
  bool                                   is_distribution_;                     // True if object is a distribution
};

#endif
