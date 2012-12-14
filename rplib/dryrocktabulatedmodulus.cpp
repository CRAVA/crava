#include "rplib/dryrocktabulatedmodulus.h"

DryRockTabulatedModulus::DryRockTabulatedModulus(const double                    k,
                                                 const double                    mu,
                                                 const double                    rho,
                                                 std::vector<double>             u)
  : DryRock()
{
  k_               = k;
  mu_              = mu;
  rho_             = rho;
  u_               = u;   // u contains correlated samples used in quantiles of (k,mu,rho)
}

DryRockTabulatedModulus::DryRockTabulatedModulus(const DryRockTabulatedModulus & rhs)
  : DryRock(rhs)
{
  k_   = rhs.k_;
  mu_  = rhs.mu_;
  rho_ = rhs.rho_;
  u_   = rhs.u_;
}

DryRock *
DryRockTabulatedModulus::Clone() const {
  return new DryRockTabulatedModulus(*this);
}

void
DryRockTabulatedModulus::ComputeElasticParams(double & k, double & mu, double & rho) const {
  k   = k_;
  mu  = mu_;
  rho = rho_;
}

DryRock *
DryRockTabulatedModulus::Evolve(const std::vector<int>             & /*delta_time*/,
                                const std::vector< const DryRock * > & /*fluid*/) const {
  return new DryRockTabulatedModulus(*this);  // Evolve when model is defined.
}
