#ifndef RPLIB_DISTRIBUTIONSFLUIDMIX_H
#define RPLIB_DISTRIBUTIONSFLUIDMIX_H

#include "rplib/distributionsfluid.h"
#include "rplib/fluidmixed.h"
#include "rplib/distributionsfluidmixevolution.h"

class DistributionsFluidMix : public DistributionsFluid {
public:

  DistributionsFluidMix(std::vector< DistributionsFluid * >            distr_fluid,
                        std::vector< NRLib::Distribution<double> * >   distr_vol_frac,
                        DistributionsFluidMixEvolution               * distr_evolution = NULL)
  : DistributionsFluid()
  {
    assert(distr_fluid.size() == distr_vol_frac.size());
    distr_fluid_      = distr_fluid;
    distr_vol_frac_   = distr_vol_frac;
    distr_evolution_  = distr_evolution;
  }

  virtual ~DistributionsFluidMix(){}

  virtual Fluid * GenerateSample() const {

    size_t n_fluids = distr_fluid_.size();
    std::vector<Fluid*> fluid(n_fluids);
    std::vector<double> volume_fraction(n_fluids);

    for(size_t i = 0; i < n_fluids; ++i) {
      fluid[i] = distr_fluid_[i]->GenerateSample();
      volume_fraction[i] = distr_vol_frac_[i]->Draw();
    }
    Fluid * fluid_mixed = new FluidMixed(fluid, volume_fraction, distr_evolution_);

    // Deep copy taken by constructor of FluidMixed, hence delete fluid here:
    for(size_t i = 0; i < n_fluids; ++i)
      delete fluid[i];

    return fluid_mixed;
  }

private:
  std::vector< DistributionsFluid * >            distr_fluid_;
  std::vector< NRLib::Distribution<double> * >   distr_vol_frac_;
  DistributionsFluidMixEvolution               * distr_evolution_;
};

#endif
