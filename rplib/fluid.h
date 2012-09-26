#ifndef RPLIB_FLUID_H
#define RPLIB_FLUID_H

#include <string>
#include <vector>

// Abstract fluid class.
class Fluid {
public:

                              Fluid() {}
  virtual                     ~Fluid() {}

  virtual Fluid *             Clone()                                                               const = 0;

  void                        GetElasticParams(double& k, double& rho) const { k = k_; rho = rho_; }

  // Fluid is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Input parameters:
  //      delta_time : the set of previous and present incremental time steps
  //      fluid : the set of previous fluid samples 
  // Recommended in implementation: assert(delta_time.size() == fluid.size() + 1);
  virtual Fluid *             Evolve(const std::vector<int>             & delta_time,
                                     const std::vector< const Fluid * > & fluid)                    const = 0;

  const std::vector<double>&  GetU()                                                                const { return u_; }


protected:
  double                      k_;
  double                      rho_;
  std::vector<double>         u_;
};

#endif
