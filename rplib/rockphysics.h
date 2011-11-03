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
  void          GetSeismicParamsForBaseFluid(double & vp, double & vs, double & rho)       const { rock_->GetSeismicParams(vp, vs, rho);} 
  void          GetSeismicParamsForFluid(const std::vector< const Fluid * >   fluid, 
                                         const std::vector< double >          saturation, 
                                         double                             & vp, 
                                         double                             & vs, 
                                         double                             & rho)         const;
  void          GetPoro(double & poro)                                                     const { rock_->GetPoro(poro);} // Crava's RMISSING if rock without porosity.
  void          GetKMineral(double & k_mineral)                                            const { rock_->GetKMineral(k_mineral);}
  void          GetMeanSeismicParamsForBaseFluid(double & mean_vp,
                                                 double & mean_vs,
                                                 double & mean_rho)                              { rock_->GetMeanSeismicParams(mean_vp, mean_vs, mean_rho);}
  void          GetVarVpForBaseFluid(double & var_vp)                                            { rock_->GetVarVp(var_vp);}
  void          GetVarVsForBaseFluid(double & var_vs)                                            { rock_->GetVarVs(var_vs);}
  void          GetVarRhoForBaseFluid(double & var_rho)                                          { rock_->GetVarRho(var_rho);}
  void          GetCrCovVpVsForBaseFluid(double & cr_cov_vp_vs)                                  { rock_->GetCrCovVpVs(cr_cov_vp_vs);}
  void          GetCrCovVpRhoForBaseFluid(double & cr_cov_vp_rho)                                { rock_->GetCrCovVpRho(cr_cov_vp_rho);}
  void          GetCrCovVsRhoForBaseFluid(double & cr_cov_vs_rho)                                { rock_->GetCrCovVsRho(cr_cov_vs_rho);}
  void          GetMeanPoroForBaseFluid(double & mean_poro)                                      { rock_->GetMeanPoro(mean_poro);}// Crava's RMISSING if rock without porosity.
  void          GetVarPoroForBaseFluid(double & var_poro)                                        { rock_->GetVarPoro(var_poro);}  // Crava's RMISSING if rock without porosity.

  void          ReSample()                                                                       { rock_->ReSample();}
  void          FluidSubstitute(const std::vector< const Fluid * >   fluid_old, 
                                const std::vector< double >          saturation_old, 
                                const double                         vp_old, 
                                const double                         vs_old, 
                                const double                         rho_old, 
                                const std::vector< const Fluid * >   fluid_new, 
                                const std::vector< double >          saturation_new, 
                                double                             & vp_new, 
                                double                             & vs_new, 
                                double                             & rho_new)              const;

private:
  Rock              *      rock_;           // Rock with base fluid. 
  FluidSubstitution *      fluid_subst_;    // Fluid substitution scheme. 
  

};

#endif
