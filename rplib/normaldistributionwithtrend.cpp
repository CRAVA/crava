
#include "nrlib/random/distribution.hpp"
#include "nrlib/random/normal.hpp"
#include "nrlib/random/uniform.hpp"
#include "nrlib/trend/trend.hpp"

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

NormalDistributionWithTrend::NormalDistributionWithTrend(const NormalDistributionWithTrend & dist)
: is_shared_(dist.is_shared_),
  use_trend_cube_(dist.use_trend_cube_)
{
  gaussian_ = dist.gaussian_->Clone();
  mean_     = dist.mean_    ->Clone();
  var_      = dist.var_     ->Clone();
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

  double u = NRLib::Random::Unif01();

  double value = GetQuantileValue(u, s1, s2);

  return value;
}

double
NormalDistributionWithTrend::GetQuantileValue(double u, double s1, double s2) const
{

  // Want sample from Y(s1, s2) ~ Normal(mu(s1, s2), var(s1,s2))
  // Generate sample from Z ~ Normal(0, 1)
  // Calculate y(s1, s2) = mu(s1, s2) + z*sqrt(var(s1,s2))

  double dummy = 0;

  double z = gaussian_->Quantile(u);

  double y = mean_->GetValue(s1, s2, dummy) + z * std::sqrt(var_->GetValue(s1, s2, dummy));

  return y;
}
