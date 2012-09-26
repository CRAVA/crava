#include "rplib/distributionsfluidtabulatedvelocity.h"

DistributionsFluidTabulatedVelocity::DistributionsFluidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * density,
                                                                         double                        corr_vp_density)
: vp_(vp),
  density_(density),
  corr_vp_density_(corr_vp_density)
{

  if(vp_->GetIsDistribution() == true || density_->GetIsDistribution() == true) {
    has_distribution_ = true;
  }
  else
    has_distribution_ = false;
}

DistributionsFluidTabulatedVelocity::~DistributionsFluidTabulatedVelocity()
{
  delete vp_;
  delete density_;
}

Fluid *
DistributionsFluidTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Dummy function. Needs to be implemented
  Fluid * fluid = NULL;

  return fluid;
}

bool
DistributionsFluidTabulatedVelocity::HasDistribution() const
{
  return(has_distribution_);
}

std::vector<bool>
DistributionsFluidTabulatedVelocity::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Fluid *
DistributionsFluidTabulatedVelocity::UpdateSample(const std::vector< double > & /*corr*/,
                                                  const Fluid                 & /*fluid*/) const {

  return NULL;
}
