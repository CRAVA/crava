#include "rplib/distributionssolidtabulatedvelocity.h"

DistributionsSolidTabulatedVelocity::DistributionsSolidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * vs,
                                                                         const DistributionWithTrend * density,
                                                                         double                        corr_vp_vs,
                                                                         double                        corr_vp_density,
                                                                         double                        corr_vs_density)
: vp_(vp),
  vs_(vs),
  density_(density),
  corr_vp_vs_(corr_vp_vs),
  corr_vp_density_(corr_vp_density),
  corr_vs_density_(corr_vs_density)
{
}

DistributionsSolidTabulatedVelocity::~DistributionsSolidTabulatedVelocity()
{
}

Solid *
DistributionsSolidTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Dummy function. Needs to be implemented
  Solid * solid = NULL;

  return solid;
}

bool
DistributionsSolidTabulatedVelocity::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsSolidTabulatedVelocity::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}

Solid *
DistributionsSolidTabulatedVelocity::UpdateSample(const std::vector< double > &/*corr*/,
                                                  const Solid                 & /*solid*/) const {

  return NULL;
}
