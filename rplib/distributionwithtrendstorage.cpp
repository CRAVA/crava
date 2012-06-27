
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"


DistributionWithTrendStorage::DistributionWithTrendStorage()
: is_gaussian_(false)
{
}

DistributionWithTrendStorage::DistributionWithTrendStorage(double value)
: is_gaussian_(false)
{
  distribution_ = new NRLib::Delta();
  mean_         = new NRLib::TrendConstantStorage(value);
  variance_     = new NRLib::TrendConstantStorage(0);
}

DistributionWithTrendStorage::DistributionWithTrendStorage(const NRLib::TrendStorage * trend)
: is_gaussian_(false)
{
  distribution_ = new NRLib::Delta();
  mean_         = trend;
  variance_     = new NRLib::TrendConstantStorage(0);
}

DistributionWithTrendStorage::DistributionWithTrendStorage(NRLib::Distribution<double>       * distribution,
                                                           const NRLib::TrendStorage         * mean,
                                                           const NRLib::TrendStorage         * variance,
                                                           bool                                is_gaussian)
: is_gaussian_(is_gaussian)
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

const NRLib::TrendStorage *
DistributionWithTrendStorage::CloneVariance()
{
  NRLib::TrendStorage * cloned_variance = variance_->Clone();
  return(cloned_variance);
}

const DistributionWithTrend *
DistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                            const std::vector<std::string>          & trend_cube_parameters,
                                                            const std::vector<std::vector<double> > & trend_cube_sampling,
                                                            std::string                             & errTxt)
{
  NRLib::Trend * mean_trend      = mean_              ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
  NRLib::Trend * variance_trend  = variance_          ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

  DistributionWithTrend * dist_with_trend             = new DistributionWithTrend(distribution_, mean_trend, variance_trend);

  //distribution_, mean_trend and variance_trend are deleted in dist_with_trend

  distribution_ = NULL;

  return(dist_with_trend);
}

const bool
DistributionWithTrendStorage::GetIsGaussian() const
{
  return(is_gaussian_);
}
