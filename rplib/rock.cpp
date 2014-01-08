#include "rplib/rock.h"

Rock::Rock()
{
  vp_ = vs_ = rho_ = 0.0;
}


Rock::~Rock()
{
}

void
Rock::GetSeismicParams(double & vp, double & vs, double & rho) const
{
  vp  = vp_;
  vs  = vs_;
  rho = rho_;
}
