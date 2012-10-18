#include "rplib/distributionsrocktabulatedmodulus.h"
#include "rplib/rocktabulatedvelocity.h"
#include "rplib/pdf3d.h"
#include "rplib/demmodelling.h"

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
  std::vector<double> u(3);

  for(int i=0; i<3; i++)
    u[i] = NRLib::Random::Unif01();

  Rock * new_rock = GetSample(u, trend_params);

  return new_rock;
}

Rock *
DistributionsRockTabulatedModulus::GetSample(const std::vector<double> & u,
                                             const std::vector<double> & trend_params) const
{
  std::vector<double> sample;

  sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double sample_bulk    = sample[0];
  double sample_shear   = sample[1];
  double sample_density = sample[2];

  double sample_vp;
  double sample_vs;

  DEMTools::CalcSeismicParamsFromElasticParams(sample_bulk, sample_shear, sample_density, sample_vp, sample_vs);

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
  bool has_distribution = false;

  if(bulk_modulus_->GetIsDistribution() == true || shear_modulus_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}

std::vector<bool>
DistributionsRockTabulatedModulus::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> bulk_trend    = bulk_modulus_ ->GetUseTrendCube();
  std::vector<bool> shear_trend   = shear_modulus_->GetUseTrendCube();
  std::vector<bool> density_trend = density_      ->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(bulk_trend[i] == true || shear_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;

}

bool
DistributionsRockTabulatedModulus::GetIsOkForBounding() const
{
  bool is_ok_for_bounding = false;

  bool              has_distribution = HasDistribution();
  std::vector<bool> has_trend        = HasTrend();

  if(has_distribution == false && has_trend[0] == false && has_trend[1] == false)
    is_ok_for_bounding = true;

  return(is_ok_for_bounding);
}

Rock *
DistributionsRockTabulatedModulus::UpdateSample(double                      corr_param,
                                                bool                        param_is_time,
                                                const std::vector<double> & trend,
                                                const Rock                * sample)       const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time, sample->GetAlpha());
  Rock * updated_sample = GetSample(u, trend);

  return updated_sample;
}
