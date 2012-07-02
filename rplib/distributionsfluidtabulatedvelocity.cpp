#include "rplib/distributionsfluidtabulatedvelocity.h"

DistributionsFluidTabulatedVelocity::DistributionsFluidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * density,
                                                                         const DistributionWithTrend * corr_vp_density)
: vp_(vp),
  density_(density),
  corr_vp_density_(corr_vp_density)
{
}

DistributionsFluidTabulatedVelocity::~DistributionsFluidTabulatedVelocity()
{
}

Fluid *
DistributionsFluidTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Dummy function. Needs to be implemented
  Fluid * fluid = NULL;

  return fluid;
}

