#include <assert.h>
#include "rplib/rockphysics.h"


RockPhysics::RockPhysics(Rock * rock, FluidSubstitution * fluid_subst) 
: rock_(rock), fluid_subst_(fluid_subst) 
{
}

RockPhysics::~RockPhysics()
{
}


void          
RockPhysics::GetSeismicParamsForFluid(const std::vector< const Fluid * >   fluid, 
                                      const std::vector< double >          saturation, 
                                      double                             & vp, 
                                      double                             & vs, 
                                      double                             & rho) const
{
  assert( fluid.size() == saturation.size() );
  double vp_base, vs_base, rho_base;
  rock_->GetSeismicParams(vp_base, vs_base, rho_base);
  std::vector< const Fluid * > fluid_base(1, rock_->GetBaseFluid());
  const std::vector< double > saturation_base(1, 1.0);
  fluid_subst_->DoSubstitute(rock_, fluid_base, saturation_base, vp_base, vs_base, rho_base, 
                                    fluid,      saturation,      vp,      vs,      rho);
}

void 
RockPhysics::FluidSubstitute(const std::vector< const Fluid * >   fluid_old, 
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
  fluid_subst_->DoSubstitute(rock_, fluid_old, saturation_old, vp_old, vs_old, rho_old, 
                                    fluid_new, saturation_new, vp_new, vs_new, rho_new);
}





