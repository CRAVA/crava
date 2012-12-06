
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"
#include "rplib/normaldistributionwithtrend.h"
#include "rplib/betadistributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"
#include "nrlib/iotools/stringtools.hpp"

DistributionWithTrendStorage::DistributionWithTrendStorage()
{
}

DistributionWithTrendStorage::~DistributionWithTrendStorage()
{
}

//--------------------------------------------------------------//

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage()
: is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(double mean,
                                                                     bool   is_shared)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
  //Use this constructor when mean_ is a double

  mean_                    = new NRLib::TrendConstantStorage(mean);
  distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                     bool                        is_shared)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
  mean_                    = mean;
  distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(const DeltaDistributionWithTrendStorage & dist)
: is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_)
{
  mean_ = dist.mean_->Clone();

 if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::~DeltaDistributionWithTrendStorage()
{
  delete mean_;

  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
DeltaDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                 const std::vector<std::string>          & trend_cube_parameters,
                                                                 const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                 std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {  //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend = mean_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

    distribution_with_trend_= new DeltaDistributionWithTrend(mean_trend, is_shared_);
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
DeltaDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

//--------------------------------------------------------------//

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage()
: is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
}

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                       const NRLib::TrendStorage * variance,
                                                                       bool                        is_shared)
: is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
  mean_                     = mean;
  variance_                 = variance;
  distribution_with_trend_  = NULL;
}

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage(const NormalDistributionWithTrendStorage & dist)
: is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_)
{
  mean_     = dist.mean_    ->Clone();
  variance_ = dist.variance_->Clone();

  if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

NormalDistributionWithTrendStorage::~NormalDistributionWithTrendStorage()
{
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
NormalDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                  const std::vector<std::string>          & trend_cube_parameters,
                                                                  const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                  std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {     //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend      = mean_    ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_trend  = variance_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

    distribution_with_trend_= new NormalDistributionWithTrend(mean_trend, variance_trend, is_shared_);
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
NormalDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}
//--------------------------------------------------------------//

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage()
: is_shared_(false),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
}

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                   const NRLib::TrendStorage * variance,
                                                                   const double              & lower_limit,
                                                                   const double              & upper_limit,
                                                                   bool                        is_shared)
: lower_limit_(lower_limit),
  upper_limit_(upper_limit),
  is_shared_(is_shared),
  vintage_year_(1),
  one_year_correlation_(1.0)
{
  mean_                     = mean;
  variance_                 = variance;
  distribution_with_trend_  = NULL;
}

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage(const BetaDistributionWithTrendStorage & dist)
: lower_limit_(dist.upper_limit_),
  upper_limit_(dist.lower_limit_),
  is_shared_(dist.is_shared_),
  vintage_year_(dist.vintage_year_),
  one_year_correlation_(dist.one_year_correlation_)
{
  mean_     = dist.mean_    ->Clone();
  variance_ = dist.variance_->Clone();

 if(dist.distribution_with_trend_ != NULL)
    distribution_with_trend_  = dist.distribution_with_trend_->Clone();
  else
    distribution_with_trend_ = NULL;
}

BetaDistributionWithTrendStorage::~BetaDistributionWithTrendStorage()
{
  delete mean_;
  delete variance_;
  distribution_with_trend_ = NULL;
}

DistributionWithTrend *
BetaDistributionWithTrendStorage::GenerateDistributionWithTrend(const std::string                       & path,
                                                                const std::vector<std::string>          & trend_cube_parameters,
                                                                const std::vector<std::vector<double> > & trend_cube_sampling,
                                                                std::string                             & errTxt)
{
  if(distribution_with_trend_ == NULL) {     //Make sure shared variables are only generated one time

    NRLib::Trend * mean_trend      = mean_    ->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);
    NRLib::Trend * variance_trend  = variance_->GenerateTrend(path,trend_cube_parameters,trend_cube_sampling,errTxt);

    CheckBetaConsistency(mean_trend, variance_trend, errTxt);

    distribution_with_trend_= new BetaDistributionWithTrend(mean_trend, variance_trend, lower_limit_, upper_limit_, is_shared_);
  }

  return(distribution_with_trend_);
}

NRLib::TrendStorage *
BetaDistributionWithTrendStorage::CloneMean() const
{
  NRLib::TrendStorage * cloned_mean = mean_->Clone();
  return(cloned_mean);
}

void
BetaDistributionWithTrendStorage::CheckBetaConsistency(NRLib::Trend * mean,
                                                       NRLib::Trend * variance,
                                                       std::string  & errTxt) const
{

  if(mean->GetMinValue() < lower_limit_ || mean->GetMaxValue() > upper_limit_)
    errTxt += "The mean values of the Beta distribution must be in the interval ["+NRLib::ToString(lower_limit_)+","+NRLib::ToString(upper_limit_)+"]\n";

  double mean_value = mean    ->GetTrendElement(0,0,0);
  double var_value  = variance->GetTrendElement(0,0,0);

  double diff = upper_limit_-lower_limit_;
  double moved_mean = (mean_value-lower_limit_) / diff;
  double moved_var  = var_value / std::pow(diff,2);

  double alpha = moved_mean * (moved_mean*(1-moved_mean)/moved_var - 1);
  double beta  = (1-moved_mean) * (moved_mean*(1-moved_mean)/moved_var - 1);

  if(alpha < 0) {
    errTxt += "The combination mean="+NRLib::ToString(mean_value)+", variance="+NRLib::ToString(var_value)+" in [a,b]=["+NRLib::ToString(lower_limit_)+","+NRLib::ToString(upper_limit_)+"]\n";
    errTxt += "  provides alpha < 0 in the Beta distribution\n";
    errTxt += "  Check that the upper and lower limits in the Beta distribution are correct\n";
  }
  if(beta < 0) {
    errTxt += "The combination mean ="+NRLib::ToString(mean_value)+", variance="+NRLib::ToString(var_value)+" in [a,b]=["+NRLib::ToString(lower_limit_)+","+NRLib::ToString(upper_limit_)+"]\n";
    errTxt += "  provides beta < 0 in the Beta distribution\n";
    errTxt += "  Check that the upper and lower limits in the Beta distribution are correct\n";
  }
  if(moved_var > moved_mean*(1-moved_mean)) {
    errTxt += "The Beta distribution demands that v < e*(1-e)\n";
    errTxt += "  where e = (mean-a)/(b-a) and v = variance/(b-a)^2`\n";
    errTxt += "  mean="+NRLib::ToString(mean_value)+", variance="+NRLib::ToString(var_value)+" in [a,b]=["+NRLib::ToString(lower_limit_)+","+NRLib::ToString(upper_limit_)+"] found\n";
  }
}
