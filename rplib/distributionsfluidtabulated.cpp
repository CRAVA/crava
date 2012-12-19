#include "rplib/distributionsfluidtabulated.h"
#include "rplib/fluidtabulatedmodulus.h"
#include "rplib/tabulated.h"

#include "rplib/fluid.h"
#include "rplib/distributionwithtrend.h"

#include "nrlib/random/distribution.hpp"

#include <cmath>

DistributionsFluidTabulated::DistributionsFluidTabulated(const DistributionWithTrend * elastic,
                                                         const DistributionWithTrend * density,
                                                         double                        corr_elastic_density,
                                                         DEMTools::TabulatedMethod     method,
                                                         std::vector<double>         & alpha)
: DistributionsFluid()
{
  if (elastic->GetIsShared() == true)
    elastic_ = elastic;
  else
    elastic_ = elastic->Clone();

  if (density->GetIsShared() == true)
    density_ = density;
  else
    density_ = density->Clone();

  corr_elastic_density_ = corr_elastic_density;
  tabulated_method_     = method;
  alpha_                = alpha;

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(2);
  elastic_variables[0] = elastic_;
  elastic_variables[1] = density_;

  NRLib::Grid2D<double> corr_matrix(2,2,0);

  for(int i=0; i<2; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_elastic_density_;
  corr_matrix(1,0) = corr_elastic_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

}

DistributionsFluidTabulated::DistributionsFluidTabulated(const DistributionsFluidTabulated & dist)
: DistributionsFluid(dist),
  corr_elastic_density_(dist.corr_elastic_density_),
  tabulated_method_(dist.tabulated_method_)
{
  if(dist.elastic_->GetIsShared() == false)
    elastic_ = dist.elastic_->Clone();
  else
    elastic_ = dist.elastic_;

  if(dist.density_->GetIsShared() == false)
    density_ = dist.density_->Clone();
  else
    density_ = dist.density_;

  // Generate tabulated_
  std::vector<const DistributionWithTrend *> elastic_variables(2);
  elastic_variables[0] = elastic_;
  elastic_variables[1] = density_;

  NRLib::Grid2D<double> corr_matrix(2,2,0);

  for(int i=0; i<2; i++)
    corr_matrix(i,i) = 1;

  corr_matrix(0,1) = corr_elastic_density_;
  corr_matrix(1,0) = corr_elastic_density_;

  tabulated_ = new Tabulated(elastic_variables, corr_matrix);

  alpha_ = dist.alpha_;
}

DistributionsFluidTabulated::~DistributionsFluidTabulated()
{
  if(elastic_->GetIsShared() == false)
    delete elastic_;
  if(density_->GetIsShared() == false)
    delete density_;

  delete tabulated_;
}

DistributionsFluid *
DistributionsFluidTabulated::Clone() const
{
  return new DistributionsFluidTabulated(*this);
}

Fluid *
DistributionsFluidTabulated::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(2);

  for(int i=0; i<2; i++)
    u[i] = NRLib::Random::Unif01();

  Fluid * fluid = GetSample(u, trend_params);

  return fluid;
}

Fluid *
DistributionsFluidTabulated::GetSample(const std::vector<double> & u, const std::vector<double> & trend_params) const
{
  std::vector<double> seismic_sample;

  seismic_sample = tabulated_->GetQuantileValues(u, trend_params[0], trend_params[1]);

  double elastic_sample = seismic_sample[0];
  double density_sample = seismic_sample[1];
  double bulk_sample;

  if(tabulated_method_ == DEMTools::Velocity)
    bulk_sample    = density_sample * (std::pow(elastic_sample,2));
  else
    bulk_sample = elastic_sample;

  Fluid * fluid = new FluidTabulatedModulus(bulk_sample, density_sample, u);

  return fluid;
}

bool
DistributionsFluidTabulated::HasDistribution() const
{
  bool has_distribution = false;

  if(elastic_->GetIsDistribution() == true || density_->GetIsDistribution() == true)
    has_distribution = true;

  return has_distribution;
}

std::vector<bool>
DistributionsFluidTabulated::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> elastic_trend = elastic_->GetUseTrendCube();
  std::vector<bool> density_trend = density_->GetUseTrendCube();

  for(int i=0; i<2; i++) {
    if(elastic_trend[i] == true || density_trend[i] == true)
      has_trend[i] = true;
  }

 return has_trend;
}

Fluid *
DistributionsFluidTabulated::UpdateSample(double                      corr_param,
                                          bool                        param_is_time,
                                          const std::vector<double> & trend,
                                          const Fluid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Fluid * updated_sample = GetSample(u, trend);

  return updated_sample;
}
