#ifndef FLUIDSUBSTITUTION_H 
#define FLUIDSUBSTITUTION_H

#include <vector>

#include "rplib/fluid.h"
#include "rplib/rock.h"


// Abstract class containing the fluid substitution scheme.
class FluidSubstitution {
public:

  FluidSubstitution();

  virtual ~FluidSubstitution();

  void DoSubstitute(const Rock                         * rock, 
                    const std::vector< const Fluid * >   fluid_old, 
                    const std::vector< double >          saturation_old, 
                    const double                         vp_old, 
                    const double                         vs_old, 
                    const double                         rho_old, 
                    const std::vector< const Fluid * >   fluid_new, 
                    const std::vector< double >          saturation_new, 
                    double                             & vp_new, 
                    double                             & vs_new, 
                    double                             & rho_new) const;

protected:
  virtual bool InputOK(const Rock                         * rock, 
                       const std::vector< const Fluid * >   fluid_old, 
                       const std::vector< double >          saturation_old, 
                       const std::vector< const Fluid * >   fluid_new, 
                       const std::vector< double >          saturation_new) const = 0;

  virtual void DoSubstituteImpl(const Rock                         * rock, 
                                const std::vector< const Fluid * >   fluid_old, 
                                const std::vector< double >          saturation_old, 
                                const double                         vp_old, 
                                const double                         vs_old, 
                                const double                         rho_old, 
                                const std::vector< const Fluid * >   fluid_new, 
                                const std::vector< double >          saturation_new, 
                                double                             & vp_new, 
                                double                             & vs_new, 
                                double                             & rho_new) const = 0;

  void DoGassmannReuss(const Rock                         * rock, 
                       const std::vector< const Fluid * >   fluid_old, 
                       const std::vector< double >          saturation_old, 
                       const double                         vp_old, 
                       const double                         vs_old, 
                       const double                         rho_old, 
                       const std::vector< const Fluid * >   fluid_new, 
                       const std::vector< double >          saturation_new, 
                       double                             & vp_new, 
                       double                             & vs_new, 
                       double                             & rho_new) const;  // Implementation of Eq. 1.15 in Avseth, Mukerji, Mavko.


};

#endif
