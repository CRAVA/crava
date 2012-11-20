#include "rplib/distributionsfluidco2.h"

#include "rplib/distributionwithtrend.h"
#include "rplib/demmodelling.h"

#include "rplib/fluidco2.h"

DistributionsFluidCO2::DistributionsFluidCO2(const DistributionWithTrend       * distr_temperature,
                                             const DistributionWithTrend       * distr_pore_pressure,
                                             std::vector<double>               & alpha)
: DistributionsFluid()
{
  distr_temperature_   = distr_temperature;
  distr_pore_pressure_ = distr_pore_pressure;
  alpha_               = alpha;
}

DistributionsFluidCO2::DistributionsFluidCO2(const DistributionsFluidCO2 & dist)
: DistributionsFluid(dist)
{
  distr_temperature_   = dist.distr_temperature_->Clone();
  distr_pore_pressure_ = dist.distr_pore_pressure_->Clone();
  alpha_               = dist.alpha_;
}

DistributionsFluidCO2::~DistributionsFluidCO2()
{
  if(distr_temperature_->GetIsShared() == false)
    delete distr_temperature_;

  if(distr_pore_pressure_->GetIsShared() == false)
    delete distr_pore_pressure_;
}

DistributionsFluid *
DistributionsFluidCO2::Clone() const
{
  return new DistributionsFluidCO2(*this);
}

Fluid *
DistributionsFluidCO2::GenerateSample(const std::vector<double> & trend_params) const
{
  std::vector<double> u(2);
  for (size_t i = 0; i < u.size(); i++)
    u[i] = NRLib::Random::Unif01();

  Fluid * fluid = GetSample(u, trend_params);

  return  fluid;
}

bool
DistributionsFluidCO2::HasDistribution() const
{
  bool has_distribution = false;

  if(distr_temperature_->GetIsDistribution() == true || distr_pore_pressure_->GetIsDistribution() == true)
    has_distribution = true;

  return(has_distribution);
}

std::vector<bool>
DistributionsFluidCO2::HasTrend() const
{
  std::vector<bool> has_trend(2, false);

  std::vector<bool> temperature_trend = distr_temperature_  ->GetUseTrendCube();
  std::vector<bool> pressure_trend    = distr_pore_pressure_->GetUseTrendCube();

  for (size_t i = 0;  i < has_trend.size(); i++) {
    if(temperature_trend[i] == true || pressure_trend[i] == true)
      has_trend[i] = true;
  }

  return has_trend;
}

Fluid *
DistributionsFluidCO2::UpdateSample(double                      corr_param,
                                    bool                        param_is_time,
                                    const std::vector<double> & trend,
                                    const Fluid               * sample) const
{
  std::vector<double> u = sample->GetU();
  DEMTools::UpdateU(u, corr_param, param_is_time);
  Fluid * updated_sample = GetSample(u, trend);

  return updated_sample;
}

Fluid *
DistributionsFluidCO2::GetSample(const std::vector<double> & u,
                                 const std::vector<double> & trend_params) const
{
  double  temperature   = distr_temperature_  ->GetQuantileValue(u[0], trend_params[0], trend_params[1]);
  double  pore_pressure = distr_pore_pressure_->GetQuantileValue(u[1], trend_params[0], trend_params[1]);

  Fluid * fluid         = new FluidCO2(temperature, pore_pressure, u);

  return  fluid;

}
