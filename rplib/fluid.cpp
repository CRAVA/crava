#include "rplib/fluid.h"

Fluid::Fluid(const double temp, const double pore_pressure)
{
  SetCommonParams(temp, pore_pressure);
}


Fluid::~Fluid() 
{
}
