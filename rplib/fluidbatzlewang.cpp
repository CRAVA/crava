#include "rplib/fluidbatzlewang.h"

#include "rplib/demmodelling.h"

#include <cassert>

FluidBatzleWang::FluidBatzleWang(double                        salinity,
                                 double                        temp,
                                 double                        pore_pressure,
                                 const std::vector<double>   & u)
: Fluid()
{
  salinity_        = salinity;
  u_               = u; //u contains independent samples used in quantiles of (salinity,temp,pore_pressure)

  ComputeElasticParams(temp, pore_pressure);
}

FluidBatzleWang::FluidBatzleWang(const FluidBatzleWang & rhs) : Fluid(rhs)
{
  salinity_ = rhs.salinity_;
}

FluidBatzleWang::~FluidBatzleWang(){}

Fluid *
FluidBatzleWang::Clone() const {
  return new FluidBatzleWang(*this);
}

void
FluidBatzleWang::ComputeElasticParams(double temp, double pore_pressure) {
  k_   = ComputeBulkModulusOfFluidBatzleWangFromTPS(temp, pore_pressure, salinity_);
  rho_ = ComputeDensityOfFluidBatzleWangFromTPS(temp, pore_pressure, salinity_);
}

Fluid *
FluidBatzleWang::Evolve(const std::vector<int>             & /*delta_time*/,
                        const std::vector< const Fluid * > & /*fluid*/) const {
  return new FluidBatzleWang(*this);
}

double
FluidBatzleWang::ComputeBulkModulusOfFluidBatzleWangFromTPS(double temp, double pore_pressure, double salinity) const
{
  return DEMTools::CalcBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity) * 0.001; /*gpa*/
}

double
FluidBatzleWang::ComputeDensityOfFluidBatzleWangFromTPS(double temp, double pore_pressure, double salinity) const
{
  return DEMTools::CalcDensityOfBrineFromTPS(temp, pore_pressure, salinity);
}

