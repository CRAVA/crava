#include "rplib/solidtabulatedmodulus.h"

SolidTabulatedModulus::SolidTabulatedModulus(const double                    k,
                                             const double                    mu,
                                             const double                    rho,
                                             std::vector<double>             u)
  : Solid()
{
  k_               = k;
  mu_              = mu;
  rho_             = rho;
  u_               = u;   // u contains correlated samples used in quantiles of (k,mu,rho)
}

SolidTabulatedModulus::SolidTabulatedModulus(const SolidTabulatedModulus & rhs)
  : Solid(rhs)
{
  k_   = rhs.k_;
  mu_  = rhs.mu_;
  rho_ = rhs.rho_;
  u_   = rhs.u_;
}

Solid *
SolidTabulatedModulus::Clone() const {
  return new SolidTabulatedModulus(*this);
}

void
SolidTabulatedModulus::ComputeElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

Solid *
SolidTabulatedModulus::Evolve(const std::vector<int>             & /*delta_time*/,
                              const std::vector< const Solid * > & /*fluid*/) const {
  return new SolidTabulatedModulus(*this);  // Evolve when model is defined.
}
