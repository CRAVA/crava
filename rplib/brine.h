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
        DistributionsBrineEvolution * distr_evolution = NULL);

  // Copy constructor
  Brine(const Brine & rhs);

  virtual ~Brine();

  // Assignment operator, not yet implemented.
  /*Brine& operator=(const Brine& rhs);*/

  virtual Fluid * Clone() const;

  virtual void ComputeElasticParams(double temp, double pore_pressure);

  virtual void GetElasticParams(double& k, double& rho) const;

  virtual Fluid * Evolve(const std::vector<int>       & /*delta_time*/,
                         const std::vector< Fluid * > & /*fluid*/) const;

private:
  double ComputeBulkModulusOfBrineFromTPS(double temp, double pore_pressure, double salinity) const;

  double ComputeDensityOfBrineFromTPS(double temp, double pore_pressure, double salinity) const;

  double salinity_;
  DistributionsBrineEvolution * distr_evolution_; // Pointer to external object.

  double k_;
  double rho_;
};

#endif
