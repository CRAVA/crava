#ifndef MINERAL_H
#define MINERAL_H

#include "rplib/solid.h"

class Mineral : public Solid {
public:

  Mineral(const double k, const double mu, const double rho)
  : Solid() {
    k_   = k;
    mu_  = mu;
    rho_ = rho;
  }

  virtual ~Mineral(){}

  virtual void ComputeElasticParams(double & k, double & mu, double & rho) const {
    k   = k_;
    mu  = mu_;
    rho = rho_;
  }

  virtual Solid * Evolve(const std::vector<int>             & /*delta_time*/,
                         const std::vector< Solid * >       & /*fluid*/,
                         const DistributionsSolidEvolution  * /*dist_fluid_evolve*/) const {
    return new Mineral(*this);
  }

private:
  double k_;
  double mu_;
  double rho_;
};

#endif
