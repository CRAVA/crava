#include "rplib/rocktabulatedvelocity.h"

RockTabulatedVelocity::RockTabulatedVelocity(double               vp,
                                             double               vs,
                                             double               density,
                                             std::vector<double>  u)
: Rock()
{
  vp_    = vp;
  vs_    = vs;
  rho_   = density;
  u_     = u;       // u contains correlated samples used in quantiles of (vp,vs,rho)
}

RockTabulatedVelocity::RockTabulatedVelocity()
: Rock()
{
}

RockTabulatedVelocity::RockTabulatedVelocity(const RockTabulatedVelocity & old_rock)
: Rock(old_rock)
{
  vp_  = old_rock.vp_;
  vs_  = old_rock.vs_;
  rho_ = old_rock.rho_;
  u_   = old_rock.u_;
}


RockTabulatedVelocity::~RockTabulatedVelocity()
{
}


Rock *
RockTabulatedVelocity::Clone() const
{
  return new RockTabulatedVelocity(*this);
}


Rock *
RockTabulatedVelocity::Evolve(const std::vector<int>         & /*delta_time*/,
                              const std::vector< Rock * >    & /*rock*/) const
{
  return new RockTabulatedVelocity(*this);
}

void
RockTabulatedVelocity::SetPorosity(double /*porosity*/)
{
}
