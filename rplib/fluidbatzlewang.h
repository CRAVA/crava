#ifndef RPLIB_FLUID_BATZLE_WANG_H
#define RPLIB_FLUID_BATZLE_WANG_H

#include "rplib/fluid.h"

class FluidBatzleWang : public Fluid {
public:

  FluidBatzleWang(double                        salinity,
                  double                        temp,
                  double                        pore_pressure,
                  const std::vector<double>   & u);

  // Copy constructor
  FluidBatzleWang(const FluidBatzleWang & rhs);

  virtual                         ~FluidBatzleWang();

                                  // Assignment operator, not yet implemented.
                                  /*FluidBatzleWang& operator=(const FluidBatzleWang& rhs);*/

  virtual Fluid *                 Clone() const;

  void                            ComputeElasticParams(double temp, double pore_pressure);

private:
  double                          ComputeBulkModulusOfFluidBatzleWangFromTPS(double temp, double pore_pressure, double salinity) const;

  double                          ComputeDensityOfFluidBatzleWangFromTPS(double temp, double pore_pressure, double salinity) const;

  double                          salinity_;

};

#endif
