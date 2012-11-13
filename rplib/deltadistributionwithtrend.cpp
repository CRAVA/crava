
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/delta.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/deltadistributionwithtrend.h"

DeltaDistributionWithTrend::DeltaDistributionWithTrend()
: is_shared_(false)
{
}

DeltaDistributionWithTrend::DeltaDistributionWithTrend(const NRLib::Trend * mean,
                                                       bool                 shared)
: mean_(mean),
  is_shared_(shared)
{

  dirac_ = new NRLib::Delta();

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(use_trend_cube_, mean_->GetTrendDimension(), mean_->GetReference());

}

DeltaDistributionWithTrend::DeltaDistributionWithTrend(const DeltaDistributionWithTrend & dist)
// Constructor to be used with Clone()
: is_shared_(false), //Set is_shared = false as the variable is cloned
  use_trend_cube_(dist.use_trend_cube_)
{
  dirac_ = dist.dirac_->Clone();
  mean_  = dist.mean_ ->Clone();
}

DeltaDistributionWithTrend::~DeltaDistributionWithTrend()
{
  delete dirac_;
  delete mean_;
}

double
DeltaDistributionWithTrend::ReSample(double s1, double s2) const
{

  // There is no random element in the Delta distribution
  double dummy = 0;

  double value = GetQuantileValue(dummy, s1, s2);

  return(value);

}

double
DeltaDistributionWithTrend::GetQuantileValue(double /*u*/, double s1, double s2) const
{

  // Want sample from Y(s1, s2) ~ Delta(mu(s1, s2)) which is one in mu(s1, s2) and zero otherwise; hence sample is equal to mean value

  double dummy = 0;

  double y = mean_->GetValue(s1, s2, dummy);

  return(y);

}
