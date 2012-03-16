#ifndef BRINE_H
#define BRINE_H

#include "rplib/fluid.h"

class Brine : public Fluid {
public:

  // Parallel classes are DistributionsBrineEvolution and DistributionsBrineT0.
  Brine(const double temp, const double pore_pressure, const double salinity)
  : Fluid(temp, pore_pressure) {
    salinity_ = salinity;
  }

  virtual ~Brine(){}

  virtual void ComputeElasticParams(double & k, double & rho) const {
    k   = ComputeBulkModulusOfBrineFromTPS(temp_, pore_pressure_, salinity_);
    rho = ComputeDensityOfBrineFromTPS(temp_, pore_pressure_, salinity_);
  }

  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< Fluid * >       & fluid,
                         const DistributionsFluidEvolution  * dist_fluid_evolve) const {
    assert(delta_time.size() == fluid.size() + 1);

    // Temporary implementation that just copies this brine:
    Fluid * new_fluid = new Brine(temp_, pore_pressure_, salinity_);
    return new_fluid;
  }

  void SetParams(const double temp, const double pore_pressure, const double salinity) {
    SetCommonParams(temp, pore_pressure);
    salinity_ = salinity;
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
