#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend()
: share_level_(None),
  current_u_(0),  //Ok since resample is true.
  resample_(true)
{
}

DistributionWithTrend::DistributionWithTrend(const int shareLevel,bool reSample)
: share_level_(shareLevel),
  current_u_(0),  //Shaky, should not be used with reSample = false, use the one below.
  resample_(reSample)
{
}

DistributionWithTrend::DistributionWithTrend(const int shareLevel,double currentU,bool reSample)
: share_level_(shareLevel),
  current_u_(currentU),
  resample_(reSample)
{
}


DistributionWithTrend::~DistributionWithTrend()
{
}

void
DistributionWithTrend::FindUseTrendCube(std::vector<bool> & use_trend_cube,
                                        int                 dim,
                                        int                 reference)
{
  if(dim == 1) {
    use_trend_cube[reference-1] = true;
  }
  else if(dim == 2) {
    use_trend_cube[0] = true;
    use_trend_cube[1] = true;
  }
}

double
DistributionWithTrend::GetCurrentSample(const std::vector<double> & trend_params)
{
  double samples;
  samples=GetQuantileValue(current_u_, trend_params[0], trend_params[1]);
  return samples;
}
