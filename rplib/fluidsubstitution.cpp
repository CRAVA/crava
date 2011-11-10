#include <assert.h>
#include <cmath>
#include "rplib/fluidsubstitution.h"


FluidSubstitution::FluidSubstitution()
{
}

FluidSubstitution::~FluidSubstitution()
{
}

bool  
FluidSubstitution::DoSubstitute(const Rock                         * rock, 
                                const std::vector< const Fluid * > & fluid_old, 
                                const std::vector< double >        & saturation_old, 
                                const double                         vp_old, 
                                const double                         vs_old, 
                                const double                         rho_old, 
                                const std::vector< const Fluid * > & fluid_new, 
                                const std::vector< double >        & saturation_new, 
                                double                             & vp_new, 
                                double                             & vs_new, 
                                double                             & rho_new) const
{
  bool subst_ok = InputOK(rock, fluid_old, saturation_old, vp_old, vs_old, rho_old, fluid_new, saturation_new);

  if (subst_ok){
    // Reuss average for bulk modules, effective fluid density, both for old fluid mix:
    double k_fluid_old   = 0.0;
    double rho_fluid_old = 0.0;
    size_t n             = fluid_old.size();
    for (size_t i = 0; i < n; ++i){
      k_fluid_old   += fluid_old[i]->GetBulkModulus() != 0.0 ? saturation_old[i] / fluid_old[i]->GetBulkModulus() : 0.0;
      rho_fluid_old += saturation_old[i] * fluid_old[i]->GetDensity();
    }
    if (k_fluid_old != 0.0)
      k_fluid_old = 1.0 / k_fluid_old;

    // Reuss average for bulk modules, effective fluid density, both for new fluid mix:
    double k_fluid_new   = 0.0;
    double rho_fluid_new = 0.0;
    n                    = fluid_new.size();
    for (size_t i = 0; i < n; ++i){
      k_fluid_new   += fluid_new[i]->GetBulkModulus() != 0.0 ? saturation_new[i] / fluid_new[i]->GetBulkModulus() : 0.0;
      rho_fluid_new += saturation_new[i] * fluid_new[i]->GetDensity();
    }
    if (k_fluid_new != 0.0)
      k_fluid_new = 1.0 / k_fluid_new;

    double c         = 4.0 / 3.0;
    double k_sat_old = rho_old * (vp_old*vp_old - c*vs_old*vs_old);
    double g_sat_old = rho_old * vs_old * vs_old;
    double k_mineral, poro;
    rock->GetKMineral(k_mineral);
    rock->GetPoro(poro);

    // Solving equation 
    // k_sat_new / (k_mineral - k_sat_new) - l2 = r1 - r2;
    double l2 = k_mineral != k_fluid_new ? k_fluid_new / (poro * (k_mineral - k_fluid_new)) : 0.0;
    double r1 = k_mineral != k_sat_old   ? k_sat_old   / (k_mineral - k_sat_old)            : 0.0;
    double r2 = k_mineral != k_fluid_old ? k_fluid_old / (poro * (k_mineral - k_fluid_old)) : 0.0;
    double b  = r1 - r2 + l2;
    double k_sat_new = b != -1.0 ? k_mineral * b / (1.0 + b)  : -1.0;

    double g_sat_new = g_sat_old;

    subst_ok = l2 != 0.0 && 
               r1 != 0.0 && 
               r2 != 0.0 &&  
               k_sat_new > 0.0 &&
               g_sat_new > 0.0;

    if (subst_ok){
      rho_new = rho_old + poro*(rho_fluid_new - rho_fluid_old); // Must be done before vp and vs.
      vp_new  = std::sqrt( (k_sat_new + c*g_sat_new) / rho_new);
      vs_new  = std::sqrt( g_sat_new / rho_new);
    }
    else{
      vp_new  = vp_old;
      vs_new  = vs_old;
      rho_new = rho_old;
    }
  }
  return subst_ok;
}

bool 
FluidSubstitution::InputOK(const Rock                         * rock, 
                           const std::vector< const Fluid * > & fluid_old, 
                           const std::vector< double >        & saturation_old,
                           const double                         vp_old, 
                           const double                         vs_old, 
                           const double                         rho_old,
                           const std::vector< const Fluid * > & fluid_new, 
                           const std::vector< double >        & saturation_new) const
{
  double poro     = 0.0;
  double kmineral = 0.0;
  if (rock != NULL){
    rock->GetPoro(poro);
    rock->GetKMineral(kmineral);
  }
  double k_sat_old = rho_old * (vp_old*vp_old - 4.0 *vs_old*vs_old / 3.0);
  bool isok = rock != NULL && 
              !rock->IsRMissing(poro) &&
              !rock->IsRMissing(kmineral) &&
              k_sat_old > 0.0 &&
              fluid_old.size() == saturation_old.size() && 
              fluid_new.size() == saturation_new.size();
  return isok;
}
