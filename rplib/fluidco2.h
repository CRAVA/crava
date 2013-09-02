#ifndef RPLIB_FLUID_FluidCO2_H
#define RPLIB_FLUID_FluidCO2_H

#include "rplib/fluid.h"



class FluidCO2 : public Fluid {
public:

  FluidCO2(double                      temp,
           double                      pore_pressure,
           const std::vector<double> & u);

  // Copy constructor
  FluidCO2(const FluidCO2 & rhs);

  virtual                     ~FluidCO2();

                              // Assignment operator, not yet implemented.
                              /*FluidCO2& operator=(const FluidCO2& rhs);*/

  virtual Fluid *             Clone()                                                               const;

  void                        ComputeElasticParams(double temp, double pressure);

};

#endif
