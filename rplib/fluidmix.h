#ifndef RPLIB_FLUIDMIX_H
#define RPLIB_FLUIDMIX_H

#include "rplib/fluid.h"
#include "rplib/demmodelling.h"

class FluidMix : public Fluid {
public:

  FluidMix(const std::vector<Fluid*>      & fluid,
           const std::vector<double>      & volume_fraction,
           const std::vector<double>      & u,
           DEMTools::MixMethod              mix_method);

  virtual ~FluidMix();

                                    // Assignment operator.
                                    FluidMix& operator=(const FluidMix& rhs);

  virtual Fluid *                   Clone() const;

  virtual Fluid *                   Evolve(const std::vector<int>             & delta_time,
                                           const std::vector< const Fluid * > & fluid) const;

  Fluid *                           GetSubFluid(size_t i) const { return fluid_[i]; }

private:
                                    //Copy constructor for getting base class variables , used by Clone:
                                    FluidMix(const FluidMix & rhs) : Fluid(rhs) {}

  std::vector<Fluid*>               fluid_;           // Owned and deleted by this class.
  std::vector<double>               volume_fraction_;
  DEMTools::MixMethod               mix_method_;

};

#endif
