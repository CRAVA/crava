#include "rplib/distributionsfluidtabulatedvelocity.h"
#include "rplib/fluidtabulatedmodulus.h"
#include "rplib/tabulated.h"

DistributionsFluidTabulatedVelocity::DistributionsFluidTabulatedVelocity(const DistributionWithTrend * vp,
                                                                         const DistributionWithTrend * density,
                                                                         double                        corr_vp_density)
: vp_(vp),
  density_(density),
  corr_vp_density_(corr_vp_density)
{

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(2);
  elastic_variables[0] = vp_;
  elastic_variables[1] = density_;

  NRLib::Grid2D<double> corr_matrix(2,2,0);

  for(int i=0; i<2; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_vp_density_;
  corr_matrix(1,0) = corr_vp_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

}

DistributionsFluidTabulatedVelocity::~DistributionsFluidTabulatedVelocity()
{
  if(vp_->GetIsShared() == false)
    delete vp_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

Fluid *
DistributionsFluidTabulatedVelocity::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(2);

  for(int i=0; i<2; i++)
    u[i] = NRLib::Random::Unif01();

  Fluid * fluid = GetSample(u, trend_params);

  return fluid;
}

Fluid *
DistributionsFluidTabulatedVelocity::GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const
{
  std::vector<double> seismic_sample;

  seismic_sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double vp_sample      = seismic_sample[0];
  double density_sample = seismic_sample[1];

  double bulk_sample  = density_sample * (std::pow(vp_sample,2));

  Fluid * fluid = new FluidTabulatedModulus(bulk_sample, density_sample, u);

  return fluid;
}

bool
DistributionsFluidTabulatedVelocity::HasDistribution() const
{
  bool has_distribution = false;

  if(vp_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}

std::vector<bool>
DistributionsFluidTabulatedVelocity::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> vp_trend      = vp_     ->GetUseTrendCube();
  std::vector<bool> density_trend = density_->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(vp_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

 return has_trend;
}

Fluid *
DistributionsFluidTabulatedVelocity::UpdateSample(const std::vector< double > & /*corr*/,
                                                  const Fluid                 & /*fluid*/) const {

  return NULL;
}
