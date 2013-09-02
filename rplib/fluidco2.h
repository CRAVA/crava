#ifndef RPLIB_FLUID_FluidCO2_H
#define RPLIB_FLUID_FluidCO2_H

#include "rplib/fluid.h"

#include <vector>



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
  void                        DebugTest();

  double                      GetCritTemp()                                                         const { return 30.9783;}
  double                      GetCritPressure()                                                     const { return 7.3772; }
  const std::vector<double>   GetPFunc()                                                                  { std::vector<double> pfunc(4); pfunc[0] = 0.000003883701221; pfunc[1] = 0.000930453898625; pfunc[2] = 0.092759127015099, pfunc[3] = 3.480623887319762; return pfunc;}
};

#endif
