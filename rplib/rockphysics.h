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
  bool          GetSeismicParamsForFluid(const std::vector< const Fluid * > & fluid, 
                                         const std::vector< double >        & saturation, 
                                         double                             & vp, 
                                         double                             & vs, 
                                         double                             & rho)         const;
  void          GetPoro(double & poro)                                                     const { rock_->GetPoro(poro);} // Rock's rmissing if rock without porosity.
  void          GetKMineral(double & k_mineral)                                            const { rock_->GetKMineral(k_mineral);}
  void          CalculateMeanSeismicParamsForBaseFluid(double & mean_vp,
                                                       double & mean_vs,
                                                       double & mean_rho)                        { rock_->CalculateMeanSeismicParams(mean_vp, mean_vs, mean_rho);}
  void          CalculateVarVpForBaseFluid(double & var_vp)                                      { rock_->CalculateVarVp(var_vp);}
  void          CalculateVarVsForBaseFluid(double & var_vs)                                      { rock_->CalculateVarVs(var_vs);}
  void          CalculateVarRhoForBaseFluid(double & var_rho)                                    { rock_->CalculateVarRho(var_rho);}
  void          CalculateCrCovVpVsForBaseFluid(double & cr_cov_vp_vs)                            { rock_->CalculateCrCovVpVs(cr_cov_vp_vs);}
  void          CalculateCrCovVpRhoForBaseFluid(double & cr_cov_vp_rho)                          { rock_->CalculateCrCovVpRho(cr_cov_vp_rho);}
  void          CalculateCrCovVsRhoForBaseFluid(double & cr_cov_vs_rho)                          { rock_->CalculateCrCovVsRho(cr_cov_vs_rho);}
  void          CalculateMeanPoroForBaseFluid(double & mean_poro)                                { rock_->CalculateMeanPoro(mean_poro);}// Rock's rmissing if rock without porosity.
  void          CalculateVarPoroForBaseFluid(double & var_poro)                                  { rock_->CalculateVarPoro(var_poro);}  // Rock's rmissing if rock without porosity.

  void          ReSample()                                                                       { rock_->ReSample();}
  bool          FluidSubstitute(const std::vector< const Fluid * > & fluid_old, 
                                const std::vector< double >        & saturation_old, 
                                const double                         vp_old, 
                                const double                         vs_old, 
                                const double                         rho_old, 
                                const std::vector< const Fluid * > & fluid_new, 
                                const std::vector< double >        & saturation_new, 
                                double                             & vp_new, 
                                double                             & vs_new, 
                                double                             & rho_new)              const;

  // Diagenesis of rock properties.
  // Return values are seismic parameters and porosity (may be Rock's rmissing) for the times specified in the input.
  void          DoDiagenesis(const std::vector< int > & time_lapse, // Measured in days.
                             std::vector< double >    & vp, 
                             std::vector< double >    & vs, 
                             std::vector< double >    & rho,
                             std::vector< double >    & poro)                              const { rock_->DoDiagenesis(time_lapse, vp, vs, rho, poro);}

private:
  Rock              *      rock_;           // Rock with base fluid. 
  FluidSubstitution *      fluid_subst_;    // Fluid substitution scheme. 
  

};

#endif
