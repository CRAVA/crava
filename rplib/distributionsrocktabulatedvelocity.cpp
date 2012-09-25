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

  /*std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = vp_;
  elastic_variables[1] = vs_;
  elastic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_elast1_elast2;
  corr_matrix(1,0) = corr_elast1_elast2;
  corr_matrix(0,2) = corr_elast1_dens;
  corr_matrix(2,0) = corr_elast1_dens;
  corr_matrix(1,2) = corr_elast2_dens;
  corr_matrix(2,1) = corr_elast2_dens;

  Tabulated distr_tab(elastic_variables, corr_matrix);*/
}

DistributionsRockTabulatedVelocity::~DistributionsRockTabulatedVelocity()
{
}

Rock *
DistributionsRockTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{

  Rock * new_rock = NULL; //new RockTabulatedVelocity();

  return new_rock;
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
