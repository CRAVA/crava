#include "rplib/distributionsrocktabulatedmodulus.h"
#include "rplib/rocktabulatedvelocity.h"
#include "rplib/pdf3d.h"

#include "nrlib/grid/grid2d.hpp"

DistributionsRockTabulatedModulus::DistributionsRockTabulatedModulus(const DistributionWithTrend * bulk_modulus,
                                                                     const DistributionWithTrend * shear_modulus,
                                                                     const DistributionWithTrend * density,
                                                                     double                        corr_bulk_shear,
                                                                     double                        corr_bulk_density,
                                                                     double                        corr_shear_density)
: bulk_modulus_(bulk_modulus),
  shear_modulus_(shear_modulus),
  density_(density),
  corr_bulk_shear_(corr_bulk_shear),
  corr_bulk_density_(corr_bulk_density),
  corr_shear_density_(corr_shear_density)
{

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> seismic_variables(3);
  seismic_variables[0] = bulk_modulus_;
  seismic_variables[1] = shear_modulus_;
  seismic_variables[2] = density_;

  NRLib::Grid2D<double> corr_matrix(3,3,0);

  for(int i=0; i<3; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_bulk_shear_;
  corr_matrix(1,0) = corr_bulk_shear_;
  corr_matrix(0,2) = corr_bulk_density_;
  corr_matrix(2,0) = corr_bulk_density_;
  corr_matrix(1,2) = corr_shear_density_;
  corr_matrix(2,1) = corr_shear_density_;

  tabulated_ = new Tabulated(seismic_variables, corr_matrix);

   // Find has_distribution_
  if(bulk_modulus_->GetIsDistribution() == true || shear_modulus_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution_ = true;
  else
    has_distribution_ = false;

  // Find has_trend_
  std::vector<bool> bulk_trend    = bulk_modulus_ ->GetUseTrendCube();
  std::vector<bool> shear_trend   = shear_modulus_->GetUseTrendCube();
  std::vector<bool> density_trend = density_      ->GetUseTrendCube();

  has_trend_.resize(2);
  for(int i=0; i<2; i++) {
    if(bulk_trend[i] == true || shear_trend[i] == true || density_trend[i] == true)
      has_trend_[i] = true;
    else
      has_trend_[i] = false;
  }

}

DistributionsRockTabulatedModulus::~DistributionsRockTabulatedModulus()
{
  if(bulk_modulus_->GetIsShared() == false)
    delete bulk_modulus_;

  if(shear_modulus_->GetIsShared() == false)
    delete shear_modulus_;

  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

Rock *
DistributionsRockTabulatedModulus::GenerateSample(const std::vector<double> & trend_params) const
{

  std::vector<double> u;
  std::vector<double> sample;

  sample = tabulated_->GenerateSample(u, trend_params[0], trend_params[1]);

  double sample_bulk    = sample[0];
  double sample_shear   = sample[1];
  double sample_density = sample[2];

  double sample_vp = std::sqrt(sample_bulk/sample_density + 4/3 * sample_shear/sample_density);
  double sample_vs = std::sqrt(sample_shear/sample_density);

  Rock * new_rock = new RockTabulatedVelocity(sample_vp, sample_vs, sample_density, u);

  return new_rock;
}

std::vector<double>
DistributionsRockTabulatedModulus::GetExpectation(const std::vector<double> & /*trend_params*/) const
{
  std::vector<double> dummy;
  return(dummy);
}

NRLib::Grid2D<double>
DistributionsRockTabulatedModulus::GetCovariance(const std::vector<double> & /*trend_params*/) const
{
  NRLib::Grid2D<double> dummy;
  return(dummy);
}

Pdf3D *
DistributionsRockTabulatedModulus::GeneratePdf() const
{
  Pdf3D * dummy = NULL;
  return(dummy);
}

bool
DistributionsRockTabulatedModulus::HasDistribution() const
{
  return(has_distribution_);
}

std::vector<bool>
DistributionsRockTabulatedModulus::HasTrend() const
{
  return(has_trend_);
}

bool
DistributionsRockTabulatedModulus::GetIsOkForBounding() const
{
  bool is_ok_for_bounding;

  if(has_distribution_ == false && has_trend_[0] == false && has_trend_[1] == false)
    is_ok_for_bounding = true;
  else
    is_ok_for_bounding = false;

  return(is_ok_for_bounding);
}
