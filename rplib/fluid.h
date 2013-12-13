#ifndef RPLIB_FLUID_H
#define RPLIB_FLUID_H

#include <string>
#include <vector>

// Abstract fluid class.
class Fluid {
public:

  Fluid();
  virtual                     ~Fluid();

  virtual Fluid            *  Clone()                                                               const = 0;

  void                        GetElasticParams(double & k,
                                               double & rho) const;

  const std::vector<double>&  GetU()                                                                const { return u_; }


protected:
  double                      k_;
  double                      rho_;
  std::vector<double>         u_;
};

#endif
