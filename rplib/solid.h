#ifndef RPLIB_SOLID_H
#define RPLIB_SOLID_H

#include <vector>

// Abstract solid class.
class Solid {
public:

                             Solid(){}
  virtual                    ~Solid() {}

  virtual Solid *            Clone()                                                  const = 0;

  // Solid is an abstract class, hence pointer must be used in Evolve.
  // Allocated memory (using new) MUST be deleted by caller.
  // Input parameters:
  //      delta_time : the set of previous and present incremental time steps
  //      solid : the set of previous fluid samples
  // Recommended in implementation: assert(delta_time.size() == solid.size() + 1);
  virtual Solid *            Evolve(const std::vector<int>             & delta_time,
                                    const std::vector< const Solid * > & solid)       const = 0;

  void                       GetElasticParams(double & k, double & mu, double & rho)  const {
                              k = k_; mu  = mu_; rho = rho_;
                             }

  const std::vector<double>& GetU()                                                   const { return u_; }

protected:
  double                      k_;
  double                      mu_;
  double                      rho_;
  std::vector<double>         u_;

};

#endif
