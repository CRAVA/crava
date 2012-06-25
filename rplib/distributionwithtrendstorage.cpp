
#include "rplib/distributionwithtrendstorage.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"


DistributionWithTrendStorage::DistributionWithTrendStorage()
{
}

DistributionWithTrendStorage::DistributionWithTrendStorage(const double value)
{
  distribution_ = new NRLib::Delta();
  mean_         = new NRLib::TrendConstantStorage(value);
  variance_     = new NRLib::TrendConstantStorage(0);
}

DistributionWithTrendStorage::DistributionWithTrendStorage(const NRLib::TrendStorage * trend)
{
  distribution_ = new NRLib::Delta();
  mean_         = trend;
  variance_     = new NRLib::TrendConstantStorage(0);
}

DistributionWithTrendStorage::DistributionWithTrendStorage(const NRLib::Distribution<double> * distribution,
                                                           const NRLib::TrendStorage         * mean,
                                                           const NRLib::TrendStorage         * variance)
{
  distribution_ = distribution;
  mean_         = mean;
  variance_     = variance;
}

DistributionWithTrendStorage::~DistributionWithTrendStorage()
{
  delete distribution_;
  delete mean_;
  delete variance_;
}

const NRLib::TrendStorage *
DistributionWithTrendStorage::CloneMean()
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}
