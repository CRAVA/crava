#include "rplib/distributionssolidtabulated.h"

DistributionsSolidTabulated::DistributionsSolidTabulated(const DistributionWithTrend * vp,
                                                         const DistributionWithTrend * vs,
                                                         const DistributionWithTrend * density,
                                                         const DistributionWithTrend * corr_vp_vs,
                                                         const DistributionWithTrend * corr_vp_density,
                                                         const DistributionWithTrend * corr_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  corr_vp_vs_(corr_vp_vs),
  corr_vp_density_(corr_vp_density),
  corr_vs_density_(corr_vs_density)
{
}

DistributionsSolidTabulated::~DistributionsSolidTabulated()
{
}

Solid *
DistributionsSolidTabulated::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Dummy function. Needs to be implemented
  Solid * solid = NULL;

  return solid;
}
