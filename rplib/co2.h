#ifndef RPLIB_CO2_H
#define RPLIB_CO2_H

#include "rplib/fluid.h"
#include "rplib/distributionsco2evolution.h"
#include "rplib/demmodelling.h"

#include <cassert>

class CO2 : public Fluid {
public:

  // Parallel classes are DistributionsCO2Evolution and DistributionsCO2.
  CO2(double   temp,
      double   pore_pressure)
  : Fluid() {
    ComputeElasticParams(temp, pore_pressure);
  }

  virtual ~CO2(){}

  virtual void ComputeElasticParams(double   temp,
                                    double   pore_pressure) {
    DEMTools::CalcCo2Prop(k_, rho_, temp, pore_pressure);
  }

  virtual void GetElasticParams(double& k, double& rho) const {
    k     = k_;
    rho   = rho_;
  }

  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Fluid * >       & fluid,
                         const DistributionsFluidEvolution  * dist_fluid_evolve) const {
    const DistributionsCO2Evolution * dist_c02_evolve = dynamic_cast<const DistributionsCO2Evolution*>(dist_fluid_evolve);
    assert(dist_c02_evolve != NULL);
    assert(delta_time.size() == fluid.size() + 1);

    // Temporary implementation that simply makes a copy, but illustrates the possible use of dist_c02_evolve:
    double example_param = dist_c02_evolve->Sample();
    example_param += 1.0; //FAKE to make avoid compiler warnings
    Fluid * new_fluid = new CO2(1.0, 1.0);
    return new_fluid;
  }

private:
  double k_;
  double rho_;
};

#endif
