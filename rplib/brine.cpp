#include "rplib/brine.h"

Brine::Brine(double                        salinity,
             double                        temp,
             double                        pore_pressure,
             DistributionsBrineEvolution * distr_evolution)
: Fluid()
{
  salinity_        = salinity;
  distr_evolution_ = distr_evolution;
  ComputeElasticParams(temp, pore_pressure);
}

Brine::Brine(const Brine & rhs) : Fluid(rhs)
{
  k_        = rhs.k_;
  rho_      = rhs.rho_;
  salinity_ = rhs.salinity_;
  distr_evolution_ = rhs.distr_evolution_;
}

Brine::~Brine(){}

Fluid *
Brine::Clone() const {
  return new Brine(*this);
}

void
Brine::ComputeElasticParams(double temp, double pore_pressure) {
  k_   = ComputeBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity_);
  rho_ = ComputeDensityOfBrineFromTPS(temp, pore_pressure, salinity_);
}

void
Brine::GetElasticParams(double& k, double& rho) const {
  k     = k_;
  rho   = rho_;
}

Fluid *
Brine::Evolve(const std::vector<int>             & /*delta_time*/,
              const std::vector< const Fluid * > & /*fluid*/) const {
  return new Brine(*this);
}

double
Brine::ComputeBulkModulusOfBrineFromTPS(double temp, double pore_pressure, double salinity) const
{
  return DEMTools::CalcBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity) * 0.001; /*gpa*/
}

double
Brine::ComputeDensityOfBrineFromTPS(double temp, double pore_pressure, double salinity) const
{
  return DEMTools::CalcDensityOfBrineFromTPS(temp, pore_pressure, salinity);
}

