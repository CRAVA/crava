#ifndef ROCKPHYSICS_H
#define ROCKPHYSICS_H

#include <vector>

#include "rplib/fluid.h"
#include "rplib/rock.h"
#include "rplib/fluidsubstitution.h"

// Rock physics class
class RockPhysics {
public:

  RockPhysics(Rock * rock, FluidSubstitution * fluid_subst);

  ~RockPhysics();

  const Fluid * GetRockBaseFluid()                                                         const { return( rock_->GetBaseFluid() );}
  void          GetSeismicParams(double & vp, double & vs, double & rho)                   const { rock_->GetSeismicParams(vp, vs, rho);} 
  void          GetSeismicParamsForFluid(const std::vector< const Fluid * >   fluid, 
                                         const std::vector< double >          saturation, 
                                         double                             & vp, 
                                         double                             & vs, 
                                         double                             & rho)         const;

  void GetPoro(double & poro)                                                              const { rock_->GetPoro(poro);} // Crava's RMISSING if rock without porosity.
  void ReSample()                                                                                { rock_->ReSample();}
  void FluidSubstitute(const std::vector< const Fluid * >   fluid_old, 
                       const std::vector< double >          saturation_old, 
                       const double                         vp_old, 
                       const double                         vs_old, 
                       const double                         rho_old, 
                       const std::vector< const Fluid * >   fluid_new, 
                       const std::vector< double >          saturation_new, 
                       double                             & vp_new, 
                       double                             & vs_new, 
                       double                             & rho_new)                        const;



  


private:
  Rock              *      rock_;           // Rock with base fluid. 
  FluidSubstitution *      fluid_subst_;    // Fluid substitution scheme. 
  

};

#endif
