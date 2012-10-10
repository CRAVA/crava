#include "rplib/fluidtabulatedmodulus.h"
#include "rplib/fluid.h"

FluidTabulatedModulus::FluidTabulatedModulus(const double                    k,
                                             const double                    rho,
                                             std::vector<double>             u)
  : Fluid()
{
  k_               = k;
  rho_             = rho;
  u_               = u;  // u contains correlated samples used in quantiles of (k,mu,rho)
}

FluidTabulatedModulus::FluidTabulatedModulus(const FluidTabulatedModulus & rhs)
  : Fluid(rhs)
{
  k_   = rhs.k_;
  rho_ = rhs.rho_;
  u_   = rhs.u_;
}

FluidTabulatedModulus::~FluidTabulatedModulus()
{
}

Fluid *
FluidTabulatedModulus::Clone() const {
  return new FluidTabulatedModulus(*this);
}

Fluid *
FluidTabulatedModulus::Evolve(const std::vector<int>             & /*delta_time*/,
                              const std::vector< const Fluid * > & /*fluid*/) const {
  return new FluidTabulatedModulus(*this);  // Evolve when model is defined.
}
