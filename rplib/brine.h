#ifndef RPLIB_BRINE_H
#define RPLIB_BRINE_H

#include "rplib/fluid.h"
#include "rplib/distributionsbrineevolution.h"
#include "rplib/demmodelling.h"

#include <cassert>

class Brine : public Fluid {
public:

  Brine(double                        salinity,
        double                        temp,
        double                        pore_pressure,
        DistributionsBrineEvolution * distr_evolution = NULL)
  : Fluid() {
    salinity_        = salinity;
    distr_evolution_ = distr_evolution;
    ComputeElasticParams(temp, pore_pressure);
  }

  // Copy constructor
  Brine(const Brine & rhs) : Fluid(rhs)
  {
    k_        = rhs.k_;
    rho_      = rhs.rho_;
    salinity_ = rhs.salinity_;
    distr_evolution_ = rhs.distr_evolution_;
  }

  virtual ~Brine(){}

  // Assignment operator, not yet implemented.
  /*Brine& operator=(const Brine& rhs);*/

  virtual Fluid * Clone() const {
    return new Brine(*this);
  }

  virtual void ComputeElasticParams(double   temp,
                                    double   pore_pressure) {
    k_   = ComputeBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity_);
    rho_ = ComputeDensityOfBrineFromTPS(temp, pore_pressure, salinity_);
  }

  virtual void GetElasticParams(double& k, double& rho) const {
    k     = k_;
    rho   = rho_;
  }

  virtual Fluid * Evolve(const std::vector<int>             & /*delta_time*/,
                         const std::vector< Fluid * >       & /*fluid*/) const {
    return new Brine(*this);
  }

private:
  double ComputeBulkModulusOfBrineFromTPS(double temp, double pore_pressure, double salinity) const {
    return DEMTools::CalcBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity) * 0.001; /*gpa*/
  }

  double ComputeDensityOfBrineFromTPS(double temp, double pore_pressure, double salinity) const {
    return DEMTools::CalcDensityOfBrineFromTPS(temp, pore_pressure, salinity);
  }

  double salinity_;
  DistributionsBrineEvolution * distr_evolution_;

  double k_;
  double rho_;
};

#endif
