#include "rplib/distributionsfluidtabulatedmodulus.h"
#include "rplib/fluidtabulatedmodulus.h"
#include "rplib/tabulated.h"
#include "rplib/demmodelling.h"

DistributionsFluidTabulatedModulus::DistributionsFluidTabulatedModulus(const DistributionWithTrend * bulk_modulus,
                                                                       const DistributionWithTrend * density,
                                                                       double                        corr_bulk_density)
: bulk_modulus_(bulk_modulus),
  density_(density),
  corr_bulk_density_(corr_bulk_density)
{

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(2);
  elastic_variables[0] = bulk_modulus_;
  elastic_variables[1] = density_;

  NRLib::Grid2D<double> corr_matrix(2,2,0);

  for(int i=0; i<2; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_bulk_density_;
  corr_matrix(1,0) = corr_bulk_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);
}

DistributionsFluidTabulatedModulus::~DistributionsFluidTabulatedModulus()
{
  if(bulk_modulus_->GetIsShared() == false)
    delete bulk_modulus_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

Fluid *
DistributionsFluidTabulatedModulus::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(2);

  for(int i=0; i<2; i++)
    u[i] = NRLib::Random::Unif01();

  Fluid * fluid = GetSample(u, trend_params);

  return fluid;
}

Fluid *
DistributionsFluidTabulatedModulus::GetSample(const std::vector<double> & u,
                                              const std::vector<double> & trend_params) const
{
  std::vector<double> sample;

  sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double bulk_sample    = sample[0];
  double density_sample = sample[1];

  Fluid * fluid = new FluidTabulatedModulus(bulk_sample, density_sample, u);

  return fluid;
}

bool
DistributionsFluidTabulatedModulus::HasDistribution() const
{
  bool has_distribution = false;

  if(bulk_modulus_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}

std::vector<bool>
DistributionsFluidTabulatedModulus::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> k_trend       = bulk_modulus_  ->GetUseTrendCube();
  std::vector<bool> density_trend = density_       ->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(k_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

Fluid *
DistributionsFluidTabulatedModulus::UpdateSample(double                      corr_param,
                                                 bool                        param_is_time,
                                                 const std::vector<double> & trend,
                                                 const Fluid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Fluid * updated_sample = GetSample(u, trend);

  return updated_sample;
}
