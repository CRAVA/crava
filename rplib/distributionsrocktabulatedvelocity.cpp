#include "rplib/distributionsrocktabulatedvelocity.h"
#include "rplib/pdf3d.h"

#include "nrlib/grid/grid2d.hpp"

DistributionsRockTabulatedVelocity::DistributionsRockTabulatedVelocity(const DistributionWithTrend * vp,
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

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(3);
  elastic_variables[0] = vp_;
  elastic_variables[1] = vs_;
  elastic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_vp_vs_;
  corr_matrix(1,0) = corr_vp_vs_;
  corr_matrix(0,2) = corr_vp_density_;
  corr_matrix(2,0) = corr_vp_density_;
  corr_matrix(1,2) = corr_vs_density_;
  corr_matrix(2,1) = corr_vs_density_;

  tabulated_ = Tabulated(elastic_variables, corr_matrix);

   // Find has_distribution_
  if(vp_->GetIsDistribution() == true || vs_->GetIsDistribution() == true || density_->GetIsDistribution() == true) {
    has_distribution_ = true;
  }
  else
    has_distribution_ = false;

  // Find has_trend_
  std::vector<bool> vp_trend      = vp_     ->GetUseTrendCube();
  std::vector<bool> vs_trend      = vs_     ->GetUseTrendCube();
  std::vector<bool> density_trend = density_->GetUseTrendCube();

  has_trend_.resize(2);
  for(int i=0; i<2; i++) {
    if(vp_trend[i] == true || vs_trend[i] == true || density_trend[i] == true)
      has_trend_[i] = true;
    else
      has_trend_[i] = false;
  }

}

DistributionsRockTabulatedVelocity::~DistributionsRockTabulatedVelocity()
{
  if(vp_->GetIsShared() == false)
    delete vp_;

  if(vs_->GetIsShared() == false)
    delete vs_;

  if(density_->GetIsShared() == false)
    delete density_;
}

Rock *
DistributionsRockTabulatedVelocity::GenerateSample(const std::vector<double> & /*trend_params*/) const
{

  /*std::vector<double> u;

  std::vector<double> sample;

  sample = tabulated_.GenerateSample(u, trend_params[0], trend_params[1]);*/

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
  return(has_distribution_);
}

std::vector<bool>
DistributionsRockTabulatedVelocity::HasTrend() const
{
  return(has_trend_);
}
