#include "rplib/solid.h"

Solid::Solid()
{
  k_   = 0;
  mu_  = 0;
  rho_ = 0;
}

Solid::~Solid()
{
}

void
Solid::GetElasticParams(double & k,
                        double & mu,
                        double & rho) const
{
  k   = k_;
  mu  = mu_;
  rho = rho_;
}
