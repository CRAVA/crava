#ifndef RPLIB_FLUID_H
#define RPLIB_FLUID_H

#include <string>
#include <vector>

// Abstract fluid class.
class Fluid {
public:

  Fluid() {}
  //Fluid(const Fluid& /*rhs*/){}
  virtual ~Fluid() {}

  // Assignment operator, not yet implemented.
  /*Fluid& operator=(const Fluid& rhs);*/

  virtual Fluid * Clone() const = 0;

  virtual void ComputeElasticParams(double   temp,
                                    double   pore_pressure) = 0;

  virtual void GetElasticParams(double& k, double& rho) const = 0;

  // Fluid is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Input parameters:
  //      delta_time : the set of previous and present incremental time steps
  //      fluid : the set of previous fluid samples 
  // Recommended in implementation: assert(delta_time.size() == fluid.size() + 1);
  virtual Fluid * Evolve(const std::vector<int>             & delta_time,
                         const std::vector< const Fluid * > & fluid) const = 0;


protected:
};

#endif
