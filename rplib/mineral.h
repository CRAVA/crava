#ifndef RPLIB_MINERAL_H
#define RPLIB_MINERAL_H

#include "rplib/solid.h"
#include "rplib/distributionsmineralevolution.h"

class Mineral : public Solid {
public:

  Mineral(const double k,
          const double mu,
          const double rho,
          DistributionsMineralEvolution * distr_evolution = NULL);

  // Copy constructor
  Mineral(const Mineral & rhs);

  virtual ~Mineral(){}

  // Assignment operator, not yet implemented.
  /*Mineral& operator=(const Mineral& rhs);*/

  virtual Solid * Clone() const;

  virtual void ComputeElasticParams(double & k, double & mu, double & rho) const;

  virtual Solid * Evolve(const std::vector<int>       & /*delta_time*/,
                         const std::vector< Solid * > & /*fluid*/) const;

private:
  DistributionsMineralEvolution * distr_evolution_;

  double k_;
  double mu_;
  double rho_;
};

#endif
