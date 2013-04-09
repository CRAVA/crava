
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"

#include "nrlib/trend/trend.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"

DeltaDistributionWithTrend::DeltaDistributionWithTrend()
: DistributionWithTrend()
{
}

DeltaDistributionWithTrend::DeltaDistributionWithTrend(const NRLib::Trend * mean,
                                                       int                  shared)
: DistributionWithTrend(shared,true)
{
  mean_ = mean->Clone();

  dirac_ = new NRLib::Delta();

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(use_trend_cube_, mean_->GetTrendDimension(), mean_->GetReference());

}

DeltaDistributionWithTrend::DeltaDistributionWithTrend(const DeltaDistributionWithTrend & dist)
  : DistributionWithTrend(dist.share_level_,dist.current_u_,dist.resample_),
  use_trend_cube_(dist.use_trend_cube_)
{
  dirac_ = dist.dirac_->Clone();
  mean_  = dist.mean_ ->Clone();
}

DeltaDistributionWithTrend::~DeltaDistributionWithTrend()
{
  delete mean_;
  delete dirac_;
}

double
DeltaDistributionWithTrend::ReSample(double s1, double s2)
{

  // There is no random element in the Delta distribution
  double dummy = 0;

  double value = GetQuantileValue(dummy, s1, s2);

  return value;

}

double
DeltaDistributionWithTrend::GetQuantileValue(double /*u*/, double s1, double s2)
{

  // Want sample from Y(s1, s2) ~ Delta(mu(s1, s2)) which is one in mu(s1, s2) and zero otherwise; hence sample is equal to mean value

  double dummy = 0;

  double y = mean_->GetValue(s1, s2, dummy);

  return(y);

}
