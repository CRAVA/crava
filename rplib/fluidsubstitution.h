#ifndef FLUIDSUBSTITUTION_H 
#define FLUIDSUBSTITUTION_H

#include <vector>

#include "rplib/fluid.h"
#include "rplib/rock.h"


// Base class containing the fluid substitution scheme.
class FluidSubstitution {
public:

  FluidSubstitution();

  virtual ~FluidSubstitution();

  // Base class uses Gassmann-Reuss fluid substitution, with an
  // implementation of section 1.3.1, Eq. 1.15, 1.16 in Avseth, Mukerji, Mavko.
  // Virtuality ensures derived classes can use another fluid substitution.
  virtual bool DoSubstitute(const Rock                         * rock, 
                            const std::vector< const Fluid * > & fluid_old, 
                            const std::vector< double >        & saturation_old, 
                            const double                         vp_old, 
                            const double                         vs_old, 
                            const double                         rho_old, 
                            const std::vector< const Fluid * > & fluid_new, 
                            const std::vector< double >        & saturation_new, 
                            double                             & vp_new, 
                            double                             & vs_new, 
                            double                             & rho_new) const;

private:
  bool InputOK(const Rock                         * rock, 
               const std::vector< const Fluid * > & fluid_old, 
               const std::vector< double >        & saturation_old,
               const double                         vp_old, 
               const double                         vs_old, 
               const double                         rho_old, 
               const std::vector< const Fluid * > & fluid_new, 
               const std::vector< double >        & saturation_new) const;
};

#endif
