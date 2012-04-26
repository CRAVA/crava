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
      DistributionsCO2Evolution * distr_evolution = NULL);

  // Copy constructor
  CO2(const CO2 & rhs);

  virtual ~CO2();

  // Assignment operator, not yet implemented.
  /*CO2& operator=(const CO2& rhs);*/

  virtual Fluid * Clone() const;

  virtual void ComputeElasticParams(double temp, double pore_pressure);

  virtual void GetElasticParams(double& k, double& rho) const;

  virtual Fluid * Evolve(const std::vector<int>       & /*delta_time*/,
                         const std::vector< Fluid * > & /*fluid*/) const;

private:
  DistributionsCO2Evolution * distr_evolution_; // Pointer to external object.

  double k_;
  double rho_;
};

#endif
