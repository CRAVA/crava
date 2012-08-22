#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend(const NRLib::Distribution<double> * distr,
                                             const NRLib::Trend                * mean,
                                             const NRLib::Trend                * sd,
                                             bool                                sheared,
                                             bool                                is_distribution)
: distribution_(distr),
  mean_(mean), 
  sd_(sd), 
  is_sheared_(sheared), 
  is_distribution_(is_distribution)
{

  use_trend_cube_.resize(2);
  for(int i=0; i<2; i++)
    use_trend_cube_[i] = false;

  FindUseTrendCube(mean_->GetTrendDimension(), mean_->GetReference());
  FindUseTrendCube(sd_  ->GetTrendDimension(), sd_  ->GetReference());

}

DistributionWithTrend::~DistributionWithTrend()
{
  delete distribution_;
  delete mean_;
  delete sd_;
}

void
DistributionWithTrend::FindUseTrendCube(int dim, int reference)
{
  if(dim == 1) {
    use_trend_cube_[reference-1] = true;
  }
  else if(dim == 2) {
    use_trend_cube_[0] = true;
    use_trend_cube_[1] = true;
  }
}
