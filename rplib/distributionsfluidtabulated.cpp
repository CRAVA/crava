#include "rplib/distributionsfluidtabulated.h"

DistributionsFluidTabulated::DistributionsFluidTabulated(const DistributionWithTrend * vp,
                                                         const DistributionWithTrend * density,
                                                         const DistributionWithTrend * corr_vp_density)
: vp_(vp),
  density_(density),
  corr_vp_density_(corr_vp_density)
{
}

DistributionsFluidTabulated::~DistributionsFluidTabulated()
{
}

Fluid *
DistributionsFluidTabulated::GenerateSample() const
{
  //Dummy function. Needs to be implemented
  Fluid * fluid = NULL;

  return fluid;
}

