#include "rplib/distributionsrocktabulatedvelocity.h"
#include "rplib/pdf3d.h"

#include "nrlib/grid/grid2d.hpp"

DistributionsRockTabulatedVelocity::DistributionsRockTabulatedVelocity(const DistributionWithTrend * vp,
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

DistributionsRockTabulatedVelocity::~DistributionsRockTabulatedVelocity()
{
}

Rock *
DistributionsRockTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{
  //Dummy function. Needs to be implemented
  Rock * rock = NULL;

  return rock;
}

std::vector<double>
DistributionsRockTabulatedVelocity::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  std::vector<double> dummy;
  return(dummy);
}

NRLib::Grid2D<double>
DistributionsRockTabulatedVelocity::GetCovariance(const std::vector<double> & /*trend_params*/) const
{
  NRLib::Grid2D<double> dummy;
  return(dummy);
}

Pdf3D *
DistributionsRockTabulatedVelocity::GeneratePdf() const
{
  Pdf3D * dummy = NULL;
  return(dummy);
}

bool
DistributionsRockTabulatedVelocity::HasDistribution() const
{
  bool dummy = false;
  return(dummy);
}

std::vector<bool>
DistributionsRockTabulatedVelocity::HasTrend() const
{
  std::vector<bool> dummy(2);

  for(int i=0; i<2; i++)
    dummy[i] = false;

  return(dummy);
}
