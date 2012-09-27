#include "rplib/mineral.h"

Mineral::Mineral(const double                    k,
                 const double                    mu,
                 const double                    rho,
                 std::vector<double>             u,
                 DistributionsMineralEvolution * distr_evolution)
  : Solid()
{
  k_               = k;
  mu_              = mu;
  rho_             = rho;
  u_               = u;
  distr_evolution_ = distr_evolution;
}

Mineral::Mineral(const Mineral & rhs)
  : Solid(rhs)
{
  k_   = rhs.k_;
  mu_  = rhs.mu_;
  rho_ = rhs.rho_;
  u_   = rhs.u_;
  distr_evolution_ = rhs.distr_evolution_;
}

Solid *
Mineral::Clone() const {
  return new Mineral(*this);
}

void
Mineral::ComputeElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

Solid *
Mineral::Evolve(const std::vector<int>             & /*delta_time*/,
                const std::vector< const Solid * > & /*fluid*/) const {
  return new Mineral(*this);  // Evolve when model is defined.
}
