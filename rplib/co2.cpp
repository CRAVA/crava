#include "rplib/co2.h"

#include "rplib/demmodelling.h"

#include <cassert>

CO2::CO2(double                      temp,
         double                      pore_pressure)
: Fluid()
{
  ComputeElasticParams(temp, pore_pressure);
}

CO2::CO2(const CO2 & rhs) : Fluid(rhs)
{
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

Fluid *
CO2::Evolve(const std::vector<int>             & /*delta_time*/,
            const std::vector< const Fluid * > & /*fluid*/) const {
  return new CO2(*this);
}


