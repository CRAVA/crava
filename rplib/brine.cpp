#include "rplib/brine.h"

#include "rplib/demmodelling.h"

#include <cassert>

Brine::Brine(double                        salinity,
             double                        temp,
             double                        pore_pressure)
: Fluid()
{
  salinity_        = salinity;
  ComputeElasticParams(temp, pore_pressure);
}

Brine::Brine(const Brine & rhs) : Fluid(rhs)
{
  salinity_ = rhs.salinity_;
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

