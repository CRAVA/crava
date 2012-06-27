#ifndef RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H
#define RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H

#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/trend/trendstorage.hpp"

class DistributionWithTrendStorage
{
public:
  DistributionWithTrendStorage();

  DistributionWithTrendStorage(double value);

  DistributionWithTrendStorage(const NRLib::TrendStorage * trend);

  DistributionWithTrendStorage(NRLib::Distribution<double>       * distribution,
                               const NRLib::TrendStorage         * mean,
                               const NRLib::TrendStorage         * variance);

  ~DistributionWithTrendStorage();

  const NRLib::TrendStorage            * CloneMean();

  const DistributionWithTrend          * GenerateDistributionWithTrend(const std::string                       & path,
                                                                       const std::vector<std::string>          & trend_cube_parameters,
                                                                       const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                       std::string                             & errTxt);

private:
  NRLib::Distribution<double>          * distribution_;
  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
};

#endif
