
#include "rplib/distributionwithtrendstorage.h"
#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"
#include "rplib/normaldistributionwithtrend.h"
#include "rplib/betadistributionwithtrend.h"

#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"
#include "nrlib/trend/trendstorage.hpp"

DistributionWithTrendStorage::DistributionWithTrendStorage()
{
}

DistributionWithTrendStorage::~DistributionWithTrendStorage()
{
}

//--------------------------------------------------------------//

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage()
: is_shared_(false),
  vintage_number_(1),
  one_year_correlation_(1.0)
{
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(double mean,
                                                                     bool   is_shared)
: is_shared_(is_shared),
  vintage_number_(1),
  one_year_correlation_(1.0)
{
  //Use this constructor when mean_ is a double

  mean_                    = new NRLib::TrendConstantStorage(mean);
  distribution_with_trend_ = NULL;
}

DeltaDistributionWithTrendStorage::DeltaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                     bool                        is_shared)
: is_shared_(is_shared),
  vintage_number_(1),
  one_year_correlation_(1.0)
{
  mean_                    = mean;
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
  vintage_number_(1),
  one_year_correlation_(1.0)
{
}

NormalDistributionWithTrendStorage::NormalDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                       const NRLib::TrendStorage * variance,
                                                                       bool                        is_shared)
: is_shared_(is_shared),
  vintage_number_(1),
  one_year_correlation_(1.0)
{
  mean_                     = mean;
  variance_                 = variance;
  distribution_with_trend_  = NULL;
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
  vintage_number_(1),
  one_year_correlation_(1.0)
{
}

BetaDistributionWithTrendStorage::BetaDistributionWithTrendStorage(const NRLib::TrendStorage * mean,
                                                                   const NRLib::TrendStorage * variance,
                                                                   bool                        is_shared)
: is_shared_(is_shared),
  vintage_number_(1),
  one_year_correlation_(1.0)
{
  mean_                     = mean;
  variance_                 = variance;
  distribution_with_trend_  = NULL;
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

    distribution_with_trend_= new BetaDistributionWithTrend(mean_trend, variance_trend, is_shared_);
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

  if(mean->GetMinValue() < 0 || mean->GetMaxValue() > 1)
    errTxt += "The mean values of the Beta distribution must be in the interval [0,1]\n";

  std::string tmpErrTxt = "";

  std::vector<int> trend_size = mean->GetTrendSize();
  for(int i=0; i<trend_size[0]; i++) {
    for(int j=0; j<trend_size[1]; j++) {
      if(tmpErrTxt == "") {
        double mean_value = mean    ->GetTrendElement(i,j,0);
        double var_value  = variance->GetTrendElement(i,j,0);

        if(var_value > mean_value*(1-mean_value))
          tmpErrTxt += "The Beta distribution demands that variance < mean*(1-mean) in each trend location\n";
      }
    }
  }

  errTxt += tmpErrTxt;
}
