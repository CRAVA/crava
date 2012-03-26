#ifndef BRINE_H
#define BRINE_H

#include "rplib/fluid.h"
#include "rplib/distributionsbrineevolution.h"

class Brine : public Fluid {
public:

  // Parallel classes are DistributionsBrineEvolution and DistributionsBrine.
  Brine(const double salinity)
  : Fluid() {
    salinity_ = salinity;
  }

  virtual ~Brine(){}

  virtual void ComputeElasticParams(const double   temp,
                                    const double   pore_pressure,
                                    double       & k,
                                    double       & rho) const {
    k   = ComputeBulkModulusOfBrineFromTPS(temp, pore_pressure, salinity_);
    rho = ComputeDensityOfBrineFromTPS(temp, pore_pressure, salinity_);
  }

  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Fluid * >       & fluid,
                         const DistributionsFluidEvolution  * dist_fluid_evolve) const {
    const DistributionsBrineEvolution * dist_brine_evolve = dynamic_cast<const DistributionsBrineEvolution*>(dist_fluid_evolve);
    assert(dist_brine_evolve != NULL);
    assert(delta_time.size() == fluid.size() + 1);

    // Temporary implementation that simply makes a copy, but illustrates the possible use of dist_brine_evolve:
    double example_param = dist_brine_evolve->Sample();
    Fluid * new_fluid = new Brine(salinity_*example_param);

    return new_fluid;
  }

private:
  double ComputeBulkModulusOfBrineFromTPS(double temp, double pore_pressure, double salinity) const {
    double k   = 4.0;  //TEMPORARY
    return k;
  }
  double ComputeDensityOfBrineFromTPS(double temp, double pore_pressure, double salinity) const {
    double rho = 3.0;  //TEMPORARY
    return rho;
  }

  double salinity_;
};

#endif
