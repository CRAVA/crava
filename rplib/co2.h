#ifndef RPLIB_CO2_H
#define RPLIB_CO2_H

#include "rplib/fluid.h"
#include "rplib/distributionsco2evolution.h"
#include "rplib/demmodelling.h"

#include <cassert>

class CO2 : public Fluid {
public:

  CO2(double                      temp,
      double                      pore_pressure,
      DistributionsCO2Evolution * distr_evolution = NULL)
  : Fluid() {
    distr_evolution_ = distr_evolution;
    ComputeElasticParams(temp, pore_pressure);
  }

  // Copy constructor
  CO2(const CO2 & rhs) : Fluid(rhs)
  {
    k_   = rhs.k_;
    rho_ = rhs.rho_;
    distr_evolution_ = rhs.distr_evolution_;
  }

  virtual ~CO2(){}

  // Assignment operator, not yet implemented.
  /*CO2& operator=(const CO2& rhs);*/

  virtual Fluid * Clone() const {
    return new CO2(*this);
  }

  virtual void ComputeElasticParams(double   temp,
                                    double   pore_pressure) {
    DEMTools::CalcCo2Prop(k_, rho_, temp, pore_pressure);
  }

  virtual void GetElasticParams(double& k, double& rho) const {
    k     = k_;
    rho   = rho_;
  }

  virtual Fluid * Evolve(const std::vector<int>             & /*delta_time*/,
                         const std::vector< Fluid * >       & /*fluid*/) const {
    return new CO2(*this);
  }

private:
  DistributionsCO2Evolution * distr_evolution_;

  double k_;
  double rho_;
};

#endif
