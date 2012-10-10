#ifndef RPLIB_BRINE_H
#define RPLIB_BRINE_H

#include "rplib/fluid.h"

class Brine : public Fluid {
public:

  Brine(double                        salinity,
        double                        temp,
        double                        pore_pressure);

  // Copy constructor
  Brine(const Brine & rhs);

  virtual                         ~Brine();

                                  // Assignment operator, not yet implemented.
                                  /*Brine& operator=(const Brine& rhs);*/

  virtual Fluid *                 Clone() const;

  void                            ComputeElasticParams(double temp, double pore_pressure);

  virtual Fluid *                 Evolve(const std::vector<int>             & /*delta_time*/,
                                         const std::vector< const Fluid * > & /*fluid*/) const;

private:
  double                          ComputeBulkModulusOfBrineFromTPS(double temp, double pore_pressure, double salinity) const;

  double                          ComputeDensityOfBrineFromTPS(double temp, double pore_pressure, double salinity) const;

  double                          salinity_;

};

#endif
