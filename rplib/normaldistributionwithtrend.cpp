
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/normal.hpp"

#include "rplib/distributionwithtrend.h"
#include "rplib/normaldistributionwithtrend.h"

NormalDistributionWithTrend::NormalDistributionWithTrend()
: is_shared_(false)
{
}

NormalDistributionWithTrend::NormalDistributionWithTrend(const NRLib::Trend * mean,
                                                         const NRLib::Trend * var,
                                                         bool                 shared)
: mean_(mean),
  var_(var),
  is_shared_(shared)
{

  gaussian_ = new NRLib::Normal();

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(use_trend_cube_, mean_->GetTrendDimension(), mean_->GetReference());
  FindUseTrendCube(use_trend_cube_, var_ ->GetTrendDimension(), var_ ->GetReference());

}

NormalDistributionWithTrend::~NormalDistributionWithTrend()
{
  delete gaussian_;
  delete mean_;
  delete var_;
}

double
NormalDistributionWithTrend::ReSample(double s1, double s2) const
{

  // Want sample from Y(s1, s2) ~ Normal(mu(s1, s2), var(s1,s2))
  // Generate sample from Z ~ Normal(0, 1)
  // Calculate y(s1, s2) = mu(s1, s2) + z*sqrt(var(s1,s2))

  double dummy = 0;

  double y = mean_->GetValue(s1, s2, dummy) + gaussian_->Draw() * std::sqrt(var_->GetValue(s1, s2, dummy));

  return y;
}
