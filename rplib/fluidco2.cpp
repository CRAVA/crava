#include "rplib/fluidco2.h"

#include "rplib/demmodelling.h"

#include <cassert>

FluidCO2::FluidCO2(double                      temp,
                   double                      pore_pressure,
                   const std::vector<double> & u)
: Fluid()
{
  u_ = u; //u contains independent samples used in quantiles of (temp,pore_pressure)
  ComputeElasticParams(temp, pore_pressure);
}

FluidCO2::FluidCO2(const FluidCO2 & rhs) : Fluid(rhs)
{
}

FluidCO2::~FluidCO2(){}


Fluid *
FluidCO2::Clone() const
{
  return new FluidCO2(*this);
}

void
FluidCO2::ComputeElasticParams(double temp, double pore_pressure)
{
  DEMTools::CalcCo2Prop(k_, rho_, temp, pore_pressure);
}

Fluid *
FluidCO2::Evolve(const std::vector<int>             & /*delta_time*/,
                 const std::vector< const Fluid * > & /*fluid*/) const
{
  return new FluidCO2(*this);
}


