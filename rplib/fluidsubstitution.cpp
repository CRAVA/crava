#include <assert.h>
#include "rplib/fluidsubstitution.h"


FluidSubstitution::FluidSubstitution()
{
}

FluidSubstitution::~FluidSubstitution()
{
}

void  FluidSubstitution::DoSubstitute(const Rock                         * rock, 
                                      const std::vector< const Fluid * >   fluid_old, 
                                      const std::vector< double >          saturation_old, 
                                      const double                         vp_old, 
                                      const double                         vs_old, 
                                      const double                         rho_old, 
                                      const std::vector< const Fluid * >   fluid_new, 
                                      const std::vector< double >          saturation_new, 
                                      double                             & vp_new, 
                                      double                             & vs_new, 
                                      double                             & rho_new) const
{
  if (InputOK(rock, fluid_old, saturation_old, fluid_new, saturation_new))
    DoSubstituteImpl(rock, fluid_old, saturation_old, vp_old, vs_old, rho_old, 
                           fluid_new, saturation_new, vp_new, vs_new, rho_new);
  //else
  //  DoSimplestPossible(...);
}


