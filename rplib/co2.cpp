#include "rplib/co2.h"

CO2::CO2(double                      temp,
         double                      pore_pressure,
         DistributionsCO2Evolution * distr_evolution)
: Fluid()
{
  distr_evolution_ = distr_evolution;
  ComputeElasticParams(temp, pore_pressure);
}

CO2::CO2(const CO2 & rhs) : Fluid(rhs)
{
  k_   = rhs.k_;
  rho_ = rhs.rho_;
  distr_evolution_ = rhs.distr_evolution_;
}

CO2::~CO2(){}


Fluid *
CO2::Clone() const {
  return new CO2(*this);
}

void
CO2::ComputeElasticParams(double temp, double pore_pressure) {
  DEMTools::CalcCo2Prop(k_, rho_, temp, pore_pressure);
}

void
CO2::GetElasticParams(double& k, double& rho) const {
  k     = k_;
  rho   = rho_;
}

Fluid *
CO2::Evolve(const std::vector<int>             & /*delta_time*/,
            const std::vector< const Fluid * > & /*fluid*/) const {
  return new CO2(*this);
}


