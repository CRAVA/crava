
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"


DistributionWithTrendStorage::DistributionWithTrendStorage()
: is_gaussian_(false),
  is_sheared_(false),
  is_distribution_(false)
{
}

DistributionWithTrendStorage::DistributionWithTrendStorage(double value,
                                                           bool   is_sheared)
: is_gaussian_(false),
  is_sheared_(is_sheared),
  is_distribution_(false)
{
  //Use this constructor when mean_ is a double
  distribution_ = new NRLib::Delta();
  mean_         = new NRLib::TrendConstantStorage(value);
  variance_     = new NRLib::TrendConstantStorage(0);
  distribution_with_trend_ = NULL;
}

DistributionWithTrendStorage::DistributionWithTrendStorage(const NRLib::TrendStorage * trend,
                                                           bool                        is_sheared)
: is_gaussian_(false),
  is_sheared_(is_sheared),
  is_distribution_(false)

{
  //Use this constructor when mean is a trend, and no distribution is given
  distribution_ = new NRLib::Delta();
  mean_         = trend;
  variance_     = new NRLib::TrendConstantStorage(0);
  distribution_with_trend_ = NULL;
}

DistributionWithTrendStorage::DistributionWithTrendStorage(NRLib::Distribution<double>       * distribution,
                                                           const NRLib::TrendStorage         * mean,
                                                           const NRLib::TrendStorage         * variance,
                                                           bool                                is_gaussian,
                                                           bool                                is_sheared)
: is_gaussian_(is_gaussian),
  is_sheared_(is_sheared),
  is_distribution_(true)
{
  //Use this constructor when a distribution is given
  distribution_ = distribution;
  mean_         = mean;
  variance_     = variance;
  distribution_with_trend_ = NULL;
}

DistributionWithTrendStorage::~DistributionWithTrendStorage()
{
  delete distribution_;
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

NRLib::TrendStorage *
DistributionWithTrendStorage::CloneMean()
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

NRLib::TrendStorage *
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
  //Make sure sheared variables are only generated one time
  if(distribution_with_trend_ == NULL) {
    NRLib::Trend * mean_trend      = mean_    ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_trend  = variance_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

    distribution_with_trend_= new DistributionWithTrend(distribution_, mean_trend, variance_trend, is_sheared_, is_distribution_); //NBNB Marit: Use variance here, but sd is used in DistributionWithTrend(). Need to fix this

    distribution_ = NULL;

    // Variables are deleted in distribution_with_trend
  }

  return(distribution_with_trend_);
}

const bool
DistributionWithTrendStorage::GetIsGaussian() const
{
  return(is_gaussian_);
}

const bool
DistributionWithTrendStorage::GetIsSheared() const
{
  return(is_sheared_);
}

const bool
DistributionWithTrendStorage::GetIsDistribution() const
{
  return(is_distribution_);
}
