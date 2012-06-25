#ifndef RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H
#define RPLIB_DISTRIBUTION_WITH_TREND_STORAGE_H

#include "nrlib/random/distribution.hpp"
#include "nrlib/trend/trendstorage.hpp"

class DistributionWithTrendStorage
{
public:
  DistributionWithTrendStorage();

  DistributionWithTrendStorage(const double value);

  DistributionWithTrendStorage(const NRLib::TrendStorage * trend);

  DistributionWithTrendStorage(const NRLib::Distribution<double> * distribution,
                               const NRLib::TrendStorage         * mean,
                               const NRLib::TrendStorage         * variance);

  ~DistributionWithTrendStorage();

  const NRLib::TrendStorage            * CloneMean();

private:
  const NRLib::Distribution<double>    * distribution_;
  const NRLib::TrendStorage            * mean_;
  const NRLib::TrendStorage            * variance_;
};

#endif
