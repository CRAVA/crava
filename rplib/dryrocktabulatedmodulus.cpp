#include "rplib/dryrocktabulatedmodulus.h"

DryRockTabulatedModulus::DryRockTabulatedModulus(const double                    k,
                                                 const double                    mu,
                                                 const double                    rho,
                                                 const double                    total_porosity,
                                                 const double                    mineral_moduli_k,
                                                 std::vector<double>             u)
  : DryRock()
{
  k_               = k;
  mu_              = mu;
  rho_             = rho;

  total_porosity_  = total_porosity;
  mineral_moduli_k_= mineral_moduli_k;

  u_               = u;   // u contains correlated samples used in quantiles of (k,mu,rho) and  uncorrelated total_porosity and mineral_moduli_k.
}

DryRockTabulatedModulus::DryRockTabulatedModulus(const DryRockTabulatedModulus & rhs)
  : DryRock(rhs)
{
  k_                = rhs.k_;
  mu_               = rhs.mu_;
  rho_              = rhs.rho_;

  total_porosity_   = rhs.total_porosity_;
  mineral_moduli_k_ = rhs.mineral_moduli_k_;

  u_                = rhs.u_;
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
