#ifndef RPLIB_CO2_H
#define RPLIB_CO2_H

#include "rplib/fluid.h"

class CO2 : public Fluid {
public:

                              CO2(double                      temp,
                                  double                      pore_pressure);

                              // Copy constructor
                              CO2(const CO2 & rhs);

  virtual                     ~CO2();

                              // Assignment operator, not yet implemented.
                              /*CO2& operator=(const CO2& rhs);*/

  virtual Fluid *             Clone()                                                               const;

  void                        ComputeElasticParams(double temp, double pore_pressure);

  virtual Fluid *             Evolve(const std::vector<int>             & /*delta_time*/,
                                     const std::vector< const Fluid * > & /*fluid*/)                const;

};

#endif
