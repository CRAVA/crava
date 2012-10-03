#include "rplib/distributionwithtrend.h"


DistributionWithTrend::DistributionWithTrend()
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
