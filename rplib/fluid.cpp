#include "rplib/fluid.h"

Fluid::Fluid()
{
  k_   = 0.0;
  rho_ = 0.0;
}


Fluid::~Fluid()
{
}

void
Fluid::GetElasticParams(double & k,
                        double & rho) const
{
  k   = k_;
  rho = rho_;
}
