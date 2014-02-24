#include "rplib/dryrock.h"

DryRock::DryRock()
{
  k_                = 0.0;
  mu_               = 0.0;
  rho_              = 0.0;
  total_porosity_   = 0.0;
  mineral_moduli_k_ = 0.0;
}


DryRock::~DryRock()
{
}

void
DryRock::GetElasticParams(double & k,
                          double & mu,
                          double & rho) const
{
  k   = k_;
  mu  = mu_;
  rho = rho_;
}
